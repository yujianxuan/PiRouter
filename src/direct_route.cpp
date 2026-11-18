#include "direct_route.h"
#include "a_star_state_grid.h"
#include "encode_coord.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <set>
#include <sstream>
#include <unordered_set>

typedef std::pair<std::vector<ExitPos>, std::vector<ExitPos>> ExitBumpPos;

namespace {

struct PatternConflict {
    bool found = false;
    int x2 = 0;
    int y2 = 0;
    int layer = 0;
    int segment_index = -1;
};

struct PendingAreaNet {
    BumpPair pair;
    std::vector<Seg> legacy_segments;
    size_t team_index = 0;
    size_t member_index = 0;
};

struct PatternCandidate {
    std::vector<Seg> segments;
    double wirelength = 0.0;
    int vias = 0;
};

PatternConflict find_first_pattern_conflict(
    const std::vector<Seg>& segments,
    const RoutingCase& rc);

bool same_point(const ResultPoint& lhs, const ResultPoint& rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.layer == rhs.layer;
}

bool point_in_grid(const ResultPoint& point, const RoutingCase& rc) {
    return point.x >= 0 && point.x < rc.grid_max_width &&
           point.y >= 0 && point.y < rc.grid_max_height &&
           point.layer >= 0 && point.layer < rc.layer_count;
}

bool queen_line(const ResultPoint& start, const ResultPoint& end) {
    if (start.layer != end.layer) {
        return start.x == end.x && start.y == end.y;
    }
    const int dx = std::abs(end.x - start.x);
    const int dy = std::abs(end.y - start.y);
    return dx == 0 || dy == 0 || dx == dy;
}

std::vector<int> candidate_key(const std::vector<Seg>& segments) {
    std::vector<int> key;
    key.reserve(1 + 6 * segments.size());
    key.push_back(static_cast<int>(segments.size()));
    for (const Seg& segment : segments) {
        key.push_back(segment.first.x);
        key.push_back(segment.first.y);
        key.push_back(segment.first.layer);
        key.push_back(segment.second.x);
        key.push_back(segment.second.y);
        key.push_back(segment.second.layer);
    }
    return key;
}

bool normalize_pattern_candidate(
    std::vector<Seg>& segments,
    const ResultPoint& start,
    const ResultPoint& end,
    const RoutingCase& rc) {
    std::vector<ResultPoint> points;
    points.reserve(segments.size() + 1);
    points.push_back(start);
    for (const Seg& segment : segments) points.push_back(segment.second);

    std::vector<ResultPoint> compact_points;
    compact_points.reserve(points.size());
    for (const ResultPoint& point : points) {
        if (compact_points.empty() || !same_point(compact_points.back(), point)) {
            compact_points.push_back(point);
        }
    }
    if (compact_points.size() < 2 || !same_point(compact_points.front(), start) ||
        !same_point(compact_points.back(), end)) {
        return false;
    }

    segments.clear();
    for (size_t i = 1; i < compact_points.size(); ++i) {
        if (!point_in_grid(compact_points[i - 1], rc) ||
            !point_in_grid(compact_points[i], rc) ||
            !queen_line(compact_points[i - 1], compact_points[i])) {
            return false;
        }
        segments.emplace_back(compact_points[i - 1], compact_points[i]);
    }

    return true;
}

double pattern_wirelength(const std::vector<Seg>& segments) {
    double length = 0.0;
    for (const Seg& segment : segments) {
        length += std::hypot(
            static_cast<double>(segment.second.x - segment.first.x),
            static_cast<double>(segment.second.y - segment.first.y));
    }
    return length;
}

int pattern_vias(const std::vector<Seg>& segments) {
    int vias = 0;
    for (const Seg& segment : segments) {
        if (segment.first.x == segment.second.x &&
            segment.first.y == segment.second.y) {
            vias += std::abs(segment.second.layer - segment.first.layer);
        }
    }
    return vias;
}

void add_path_occupancy_keys(
    const std::vector<Seg>& segments,
    std::unordered_set<uint64_t>& keys) {
    for (const Seg& segment : segments) {
        const ResultPoint& start = segment.first;
        const ResultPoint& end = segment.second;
        if (start.layer != end.layer) {
            if (start.x != end.x || start.y != end.y) continue;
            for (int layer = std::min(start.layer, end.layer);
                 layer <= std::max(start.layer, end.layer); ++layer) {
                keys.insert(encode_coord(start.x, start.y, layer));
            }
            continue;
        }
        const int dx = end.x - start.x;
        const int dy = end.y - start.y;
        const int step_x = (dx > 0) - (dx < 0);
        const int step_y = (dy > 0) - (dy < 0);
        const int steps = std::max(std::abs(dx), std::abs(dy));
        int x = start.x;
        int y = start.y;
        for (int step = 0; step <= steps; ++step) {
            keys.insert(encode_coord(x, y, start.layer));
            if (step < steps && step_x != 0 && step_y != 0) {
                keys.insert(encode_coord(
                    static_cast<float>(x) + 0.5f * step_x,
                    static_cast<float>(y) + 0.5f * step_y,
                    start.layer));
            }
            x += step_x;
            y += step_y;
        }
    }
}

std::unordered_set<uint64_t> collect_net_escape_cells(
    const BumpPair& pair,
    const DRResultDict& escape_paths) {
    std::unordered_set<uint64_t> cells;
    for (const Bump* bump : {&pair.first, &pair.second}) {
        auto path_it = escape_paths.find(bump->exit_terminal.bump_name);
        if (path_it != escape_paths.end()) {
            add_path_occupancy_keys(path_it->second, cells);
        }
    }
    return cells;
}

void add_pattern_candidate(
    std::vector<PatternCandidate>& candidates,
    std::set<std::vector<int>>& seen,
    std::vector<Seg> segments,
    const ResultPoint& start,
    const ResultPoint& end,
    const RoutingCase& rc) {
    if (!normalize_pattern_candidate(segments, start, end, rc)) return;
    std::vector<int> key = candidate_key(segments);
    if (!seen.insert(std::move(key)).second) return;
    candidates.push_back({std::move(segments), 0.0, 0});
    candidates.back().wirelength = pattern_wirelength(candidates.back().segments);
    candidates.back().vias = pattern_vias(candidates.back().segments);
}

std::vector<int> offset_coordinates(
    const std::vector<int>& centers,
    int limit) {
    static constexpr int OFFSETS[] = {
        0, 2, -2, 4, -4, 8, -8, 16, -16, 32, -32,
        64, -64, 96, -96, 128, -128, 192, -192,
        256, -256, 384, -384, 512, -512,
    };
    std::set<int> unique;
    for (int center : centers) {
        for (int offset : OFFSETS) {
            const int coordinate = center + offset;
            if (coordinate >= 0 && coordinate < limit) unique.insert(coordinate);
        }
    }
    return std::vector<int>(unique.begin(), unique.end());
}

std::vector<PatternCandidate> generate_pattern_candidates(
    const PendingAreaNet& pending,
    const RoutingCase& rc,
    bool allow_layer_change) {
    const ExitTerminal& first_exit = pending.pair.first.exit_terminal;
    const ExitTerminal& second_exit = pending.pair.second.exit_terminal;
    const ResultPoint start{
        first_exit.grid_coord_x, first_exit.grid_coord_y, first_exit.layer};
    const ResultPoint end{
        second_exit.grid_coord_x, second_exit.grid_coord_y, second_exit.layer};
    std::vector<PatternCandidate> candidates;
    if (start.layer != end.layer) return candidates;

    std::set<std::vector<int>> seen;
    if (!pending.legacy_segments.empty()) {
        add_pattern_candidate(
            candidates, seen, pending.legacy_segments, start, end, rc);
    }

    PatternConflict legacy_conflict;
    if (!pending.legacy_segments.empty()) {
        legacy_conflict = find_first_pattern_conflict(
            pending.legacy_segments, rc);
    }
    std::vector<int> x_centers{start.x, end.x, (start.x + end.x) / 2};
    std::vector<int> y_centers{start.y, end.y, (start.y + end.y) / 2};
    if (legacy_conflict.found) {
        x_centers.push_back(legacy_conflict.x2 / 2);
        y_centers.push_back(legacy_conflict.y2 / 2);
    }
    for (const Seg& segment : pending.legacy_segments) {
        if (segment.first.x == segment.second.x) {
            x_centers.push_back(segment.first.x);
        }
        if (segment.first.y == segment.second.y) {
            y_centers.push_back(segment.first.y);
        }
    }

    std::vector<int> routing_layers{start.layer};
    if (allow_layer_change) {
        std::vector<int> alternatives;
        for (int layer = 0; layer < rc.layer_count; ++layer) {
            if (layer != start.layer) alternatives.push_back(layer);
        }
        std::sort(alternatives.begin(), alternatives.end(), [&](int lhs, int rhs) {
            const int lhs_vias = std::abs(lhs - start.layer) +
                                 std::abs(lhs - end.layer);
            const int rhs_vias = std::abs(rhs - start.layer) +
                                 std::abs(rhs - end.layer);
            if (lhs_vias != rhs_vias) return lhs_vias < rhs_vias;
            return lhs < rhs;
        });
        routing_layers.insert(
            routing_layers.end(), alternatives.begin(), alternatives.end());
    }

    const std::vector<int> horizontal_spines =
        offset_coordinates(y_centers, rc.grid_max_height);
    const std::vector<int> vertical_spines =
        offset_coordinates(x_centers, rc.grid_max_width);
    const int x_direction = (end.x > start.x) - (end.x < start.x);
    const int y_direction = (end.y > start.y) - (end.y < start.y);
    for (int routing_layer : routing_layers) {
        if (routing_layer != start.layer && !pending.legacy_segments.empty()) {
            std::vector<ResultPoint> legacy_points;
            std::vector<int> cumulative_steps;
            legacy_points.push_back(start);
            cumulative_steps.push_back(0);
            int total_legacy_steps = 0;
            for (const Seg& segment : pending.legacy_segments) {
                total_legacy_steps += std::max(
                    std::abs(segment.second.x - segment.first.x),
                    std::abs(segment.second.y - segment.first.y));
                legacy_points.push_back(segment.second);
                cumulative_steps.push_back(total_legacy_steps);
            }
            auto point_at_step = [&](int target_step) {
                for (size_t i = 1; i < legacy_points.size(); ++i) {
                    if (target_step > cumulative_steps[i]) continue;
                    const int local_step = target_step - cumulative_steps[i - 1];
                    const ResultPoint& segment_start = legacy_points[i - 1];
                    const ResultPoint& segment_end = legacy_points[i];
                    const int step_x = (segment_end.x > segment_start.x) -
                                       (segment_end.x < segment_start.x);
                    const int step_y = (segment_end.y > segment_start.y) -
                                       (segment_end.y < segment_start.y);
                    return ResultPoint{
                        segment_start.x + step_x * local_step,
                        segment_start.y + step_y * local_step,
                        start.layer};
                }
                return end;
            };
            static constexpr int ACCESS_STEPS[] = {
                1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
            };
            for (int start_steps : ACCESS_STEPS) {
                for (int end_steps : ACCESS_STEPS) {
                    if (start_steps + end_steps >= total_legacy_steps) continue;
                    const int end_target_step = total_legacy_steps - end_steps;
                    const ResultPoint start_portal = point_at_step(start_steps);
                    const ResultPoint end_portal = point_at_step(end_target_step);
                    const ResultPoint start_portal_on_layer{
                        start_portal.x, start_portal.y, routing_layer};
                    const ResultPoint end_portal_on_layer{
                        end_portal.x, end_portal.y, routing_layer};
                    std::vector<Seg> lifted;
                    lifted.emplace_back(start, start_portal);
                    lifted.emplace_back(start_portal, start_portal_on_layer);
                    ResultPoint previous = start_portal_on_layer;
                    for (size_t i = 1; i + 1 < legacy_points.size(); ++i) {
                        if (cumulative_steps[i] <= start_steps ||
                            cumulative_steps[i] >= end_target_step) {
                            continue;
                        }
                        const ResultPoint bend{
                            legacy_points[i].x,
                            legacy_points[i].y,
                            routing_layer};
                        lifted.emplace_back(previous, bend);
                        previous = bend;
                    }
                    lifted.emplace_back(previous, end_portal_on_layer);
                    lifted.emplace_back(end_portal_on_layer, end_portal);
                    lifted.emplace_back(end_portal, end);
                    add_pattern_candidate(
                        candidates, seen, std::move(lifted), start, end, rc);
                }
            }
        }
        const ResultPoint planar_start{start.x, start.y, routing_layer};
        const ResultPoint planar_end{end.x, end.y, routing_layer};
        auto wrap_planar = [&](std::vector<Seg> planar_segments) {
            if (routing_layer != start.layer && planar_segments.size() == 3) {
                const ResultPoint first_bend_on_layer = planar_segments[0].second;
                const ResultPoint second_bend_on_layer = planar_segments[1].second;
                const ResultPoint first_bend{
                    first_bend_on_layer.x, first_bend_on_layer.y, start.layer};
                const ResultPoint second_bend{
                    second_bend_on_layer.x, second_bend_on_layer.y, start.layer};
                return std::vector<Seg>{
                    {start, first_bend},
                    {first_bend, first_bend_on_layer},
                    {first_bend_on_layer, second_bend_on_layer},
                    {second_bend_on_layer, second_bend},
                    {second_bend, end},
                };
            }
            std::vector<Seg> full;
            full.reserve(planar_segments.size() + 2);
            if (!same_point(start, planar_start)) {
                full.emplace_back(start, planar_start);
            }
            full.insert(full.end(), planar_segments.begin(), planar_segments.end());
            if (!same_point(planar_end, end)) {
                full.emplace_back(planar_end, end);
            }
            return full;
        };

        for (int spine_y : horizontal_spines) {
            const int start_delta = std::abs(spine_y - start.y);
            const int end_delta = std::abs(spine_y - end.y);
            const int preferred_direction = x_direction == 0 ? 1 : x_direction;
            const int sign_pairs[2][2] = {
                {preferred_direction, -preferred_direction},
                {-preferred_direction, preferred_direction},
            };
            for (const auto& signs : sign_pairs) {
                const int start_sign = signs[0];
                const int end_sign = signs[1];
                const ResultPoint first_bend{
                    start.x + start_sign * start_delta, spine_y, routing_layer};
                const ResultPoint second_bend{
                    end.x + end_sign * end_delta, spine_y, routing_layer};
                if ((x_direction > 0 && first_bend.x > second_bend.x) ||
                    (x_direction < 0 && first_bend.x < second_bend.x)) {
                    continue;
                }
                add_pattern_candidate(
                    candidates, seen,
                    wrap_planar({{planar_start, first_bend},
                                 {first_bend, second_bend},
                                 {second_bend, planar_end}}),
                    start, end, rc);
            }
        }

        for (int spine_x : vertical_spines) {
            const int start_delta = std::abs(spine_x - start.x);
            const int end_delta = std::abs(spine_x - end.x);
            const int preferred_direction = y_direction == 0 ? 1 : y_direction;
            const int sign_pairs[2][2] = {
                {preferred_direction, -preferred_direction},
                {-preferred_direction, preferred_direction},
            };
            for (const auto& signs : sign_pairs) {
                const int start_sign = signs[0];
                const int end_sign = signs[1];
                const ResultPoint first_bend{
                    spine_x, start.y + start_sign * start_delta, routing_layer};
                const ResultPoint second_bend{
                    spine_x, end.y + end_sign * end_delta, routing_layer};
                if ((y_direction > 0 && first_bend.y > second_bend.y) ||
                    (y_direction < 0 && first_bend.y < second_bend.y)) {
                    continue;
                }
                add_pattern_candidate(
                    candidates, seen,
                    wrap_planar({{planar_start, first_bend},
                                 {first_bend, second_bend},
                                 {second_bend, planar_end}}),
                    start, end, rc);
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.vias != rhs.vias) return lhs.vias < rhs.vias;
        if (std::abs(lhs.wirelength - rhs.wirelength) > 1e-9) {
            return lhs.wirelength < rhs.wirelength;
        }
        if (lhs.segments.size() != rhs.segments.size()) {
            return lhs.segments.size() < rhs.segments.size();
        }
        return candidate_key(lhs.segments) < candidate_key(rhs.segments);
    });
    return candidates;
}

json result_point_json(const ResultPoint& point) {
    return json{{"x", point.x}, {"y", point.y}, {"layer", point.layer}};
}

json exit_terminal_json(const ExitTerminal& terminal) {
    return json{
        {"bump_name", terminal.bump_name},
        {"x", terminal.grid_coord_x},
        {"y", terminal.grid_coord_y},
        {"layer", terminal.layer},
    };
}

json segments_json(const std::vector<Seg>& segments) {
    json result = json::array();
    for (const Seg& segment : segments) {
        result.push_back(json{
            {"start", result_point_json(segment.first)},
            {"end", result_point_json(segment.second)},
        });
    }
    return result;
}

bool endpoint_overlap_allowed(
    int segment_index,
    int segment_count,
    uint64_t key,
    const ResultPoint& start,
    const ResultPoint& end) {
    return (segment_index == 0 && key == encode_coord(start.x, start.y, start.layer)) ||
           (segment_index == segment_count - 1 &&
            key == encode_coord(end.x, end.y, end.layer));
}

PatternConflict find_first_pattern_conflict(
    const std::vector<Seg>& segments,
    const RoutingCase& rc) {
    PatternConflict result;
    const int segment_count = static_cast<int>(segments.size());
    for (int segment_index = 0; segment_index < segment_count; ++segment_index) {
        const Seg& segment = segments[segment_index];
        const ResultPoint& start = segment.first;
        const ResultPoint& end = segment.second;
        if (start.layer != end.layer) {
            if (start.x != end.x || start.y != end.y) {
                return {true, 2 * start.x, 2 * start.y,
                        start.layer, segment_index};
            }
            const int first_layer = std::min(start.layer, end.layer);
            const int last_layer = std::max(start.layer, end.layer);
            for (int layer = first_layer; layer <= last_layer; ++layer) {
                const uint64_t key = encode_coord(start.x, start.y, layer);
                if (rc.occupied_set.count(key) &&
                    !endpoint_overlap_allowed(
                        segment_index, segment_count, key, start, end)) {
                    return {true, 2 * start.x, 2 * start.y,
                            layer, segment_index};
                }
            }
            continue;
        }
        const int dx = end.x - start.x;
        const int dy = end.y - start.y;
        const int step_x = (dx > 0) - (dx < 0);
        const int step_y = (dy > 0) - (dy < 0);
        const int total_steps = std::max(std::abs(dx), std::abs(dy));
        int x = start.x;
        int y = start.y;
        for (int step = 0; step <= total_steps; ++step) {
            const uint64_t key = encode_coord(x, y, start.layer);
            if (rc.occupied_set.count(key) &&
                !endpoint_overlap_allowed(
                    segment_index, segment_count, key, start, end)) {
                result = {true, 2 * x, 2 * y, start.layer, segment_index};
                return result;
            }
            if (step < total_steps && step_x != 0 && step_y != 0) {
                const float mid_x = static_cast<float>(x) + 0.5f * step_x;
                const float mid_y = static_cast<float>(y) + 0.5f * step_y;
                const uint64_t mid_key = encode_coord(mid_x, mid_y, start.layer);
                if (rc.occupied_set.count(mid_key)) {
                    result = {
                        true,
                        2 * x + step_x,
                        2 * y + step_y,
                        start.layer,
                        segment_index,
                    };
                    return result;
                }
            }
            x += step_x;
            y += step_y;
        }
    }
    return result;
}

class AreaPatternTracer {
public:
    AreaPatternTracer() {
        const char* output_env = std::getenv("PIROUTER_AREA_TRACE_PATH");
        const char* nets_env = std::getenv("PIROUTER_AREA_TRACE_NETS");
        if (output_env == nullptr || nets_env == nullptr) return;
        output_path_ = output_env;
        std::ifstream input(nets_env);
        std::string line;
        while (std::getline(input, line)) {
            if (!line.empty() && line[0] != '#') targets_.insert(line);
        }
        enabled_ = !targets_.empty();
        const char* stop_env = std::getenv("PIROUTER_AREA_TRACE_STOP_AFTER_PATTERN");
        stop_after_pattern_ = stop_env != nullptr && std::string(stop_env) != "0";
        root_["format"] = "PiRouter Area Pattern trace v1";
        root_["coordinate_note"] =
            "local occupied points use doubled coordinates x2/y2; overview density counts integer nodes";
        root_["window_radius"] = local_radius_;
        root_["targets"] = json::array();
        for (const std::string& target : targets_) root_["targets"].push_back(target);
        root_["nets"] = json::object();
    }

    bool enabled() const { return enabled_; }

    bool wants(const std::string& name) const { return is_target(name); }

    void capture_area_start(const std::vector<PairTeam>& teams, const RoutingCase& rc) {
        if (!enabled_) return;
        root_["grid"] = json{
            {"width", rc.grid_max_width},
            {"height", rc.grid_max_height},
            {"layers", rc.layer_count},
        };
        for (size_t team_index = 0; team_index < teams.size(); ++team_index) {
            const PairTeam& team = teams[team_index];
            for (size_t member_index = 0; member_index < team.member.size(); ++member_index) {
                const BumpPair& pair = team.member[member_index];
                const std::string& name = pair.first.net_name;
                if (!is_target(name)) continue;
                json& record = root_["nets"][name];
                record["team_index"] = team_index;
                record["source_team_id"] = team.source_team_id;
                record["member_index"] = member_index;
                record["bump_start"] = json{
                    {"x", pair.first.grid_coord_x},
                    {"y", pair.first.grid_coord_y},
                    {"layer", pair.first.layer},
                };
                record["bump_end"] = json{
                    {"x", pair.second.grid_coord_x},
                    {"y", pair.second.grid_coord_y},
                    {"layer", pair.second.layer},
                };
                record["exit_start"] = exit_terminal_json(pair.first.exit_terminal);
                record["exit_end"] = exit_terminal_json(pair.second.exit_terminal);
                record["team_bounds"] = json{
                    {"start", json::array({team.start_left_bound_x, team.start_right_bound_x,
                                             team.start_up_bound_y, team.start_down_bound_y})},
                    {"end", json::array({team.end_left_bound_x, team.end_right_bound_x,
                                           team.end_up_bound_y, team.end_down_bound_y})},
                };
                record["area_start"] = capture_snapshot(pair, rc);
            }
        }
    }

    void record_attempt(
        const BumpPair& pair,
        size_t team_index,
        size_t member_index,
        const std::string& pattern_type,
        const std::vector<Seg>& segments,
        const std::string& outcome,
        const PatternConflict& conflict,
        const RoutingCase& rc) {
        if (!is_target(pair.first.net_name)) return;
        json& record = root_["nets"][pair.first.net_name];
        record["team_index"] = team_index;
        record["member_index"] = member_index;
        record["pattern_type"] = pattern_type;
        record["pattern_segments"] = segments_json(segments);
        record["pattern_outcome"] = outcome;
        record["before_attempt"] = capture_snapshot(pair, rc);
        if (conflict.found) {
            record["first_conflict"] = json{
                {"x2", conflict.x2},
                {"y2", conflict.y2},
                {"x", conflict.x2 / 2.0},
                {"y", conflict.y2 / 2.0},
                {"layer", conflict.layer},
                {"segment_index", conflict.segment_index},
                {"local_occupied", capture_local_x2(
                    conflict.x2, conflict.y2, conflict.layer, rc)},
            };
        }
    }

    void finish_pattern_stage(const std::vector<PairTeam>& teams, const RoutingCase& rc) {
        if (!enabled_) return;
        for (const PairTeam& team : teams) {
            for (const BumpPair& pair : team.member) {
                const std::string& name = pair.first.net_name;
                if (!is_target(name)) continue;
                json& record = root_["nets"][name];
                record["after_all_patterns"] = capture_snapshot(pair, rc);
                record["astar_start_neighbors"] = astar_neighbor_snapshot(pair, rc);
            }
        }
        root_["missing_targets"] = json::array();
        for (const std::string& target : targets_) {
            if (!root_["nets"].contains(target) ||
                !root_["nets"][target].contains("pattern_outcome")) {
                root_["missing_targets"].push_back(target);
            }
        }
        std::ofstream output(output_path_);
        output << root_.dump(2) << '\n';
        output.close();
        std::cout << "[AreaTrace] wrote " << output_path_ << " for "
                  << targets_.size() << " targets" << std::endl;
        if (stop_after_pattern_) {
            std::cout << "[AreaTrace] stopping after Pattern stage as requested" << std::endl;
            std::exit(0);
        }
    }

private:
    bool enabled_ = false;
    bool stop_after_pattern_ = false;
    int local_radius_ = 60;
    std::string output_path_;
    std::set<std::string> targets_;
    json root_;

    bool is_target(const std::string& name) const {
        return enabled_ && targets_.count(name) != 0;
    }

    json capture_local_x2(int center_x2, int center_y2, int layer, const RoutingCase& rc) const {
        json points = json::array();
        const int radius_x2 = 2 * local_radius_;
        for (int x2 = center_x2 - radius_x2; x2 <= center_x2 + radius_x2; ++x2) {
            if (x2 < 0 || x2 >= 2 * rc.grid_max_width) continue;
            for (int y2 = center_y2 - radius_x2; y2 <= center_y2 + radius_x2; ++y2) {
                if (y2 < 0 || y2 >= 2 * rc.grid_max_height) continue;
                const float x = 0.5f * static_cast<float>(x2);
                const float y = 0.5f * static_cast<float>(y2);
                if (rc.occupied_set.count(encode_coord(x, y, layer))) {
                    points.push_back(json::array({x2, y2}));
                }
            }
        }
        return json{
            {"center_x2", center_x2},
            {"center_y2", center_y2},
            {"radius", local_radius_},
            {"layer", layer},
            {"points_x2", std::move(points)},
        };
    }

    json capture_density(const BumpPair& pair, const RoutingCase& rc) const {
        const int margin = 80;
        const int xmin = std::max(
            0,
            std::min(pair.first.exit_terminal.grid_coord_x,
                     pair.second.exit_terminal.grid_coord_x) - margin);
        const int xmax = std::min(
            rc.grid_max_width - 1,
            std::max(pair.first.exit_terminal.grid_coord_x,
                     pair.second.exit_terminal.grid_coord_x) + margin);
        const int ymin = std::max(
            0,
            std::min(pair.first.exit_terminal.grid_coord_y,
                     pair.second.exit_terminal.grid_coord_y) - margin);
        const int ymax = std::min(
            rc.grid_max_height - 1,
            std::max(pair.first.exit_terminal.grid_coord_y,
                     pair.second.exit_terminal.grid_coord_y) + margin);
        const int nx = 48;
        const int ny = 32;
        std::vector<int> counts(static_cast<size_t>(nx * ny), 0);
        std::vector<int> samples(static_cast<size_t>(nx * ny), 0);
        const int layer = pair.first.exit_terminal.layer;
        const int width = std::max(1, xmax - xmin + 1);
        const int height = std::max(1, ymax - ymin + 1);
        for (int x = xmin; x <= xmax; ++x) {
            const int bx = std::min(nx - 1, (x - xmin) * nx / width);
            for (int y = ymin; y <= ymax; ++y) {
                const int by = std::min(ny - 1, (y - ymin) * ny / height);
                const size_t index = static_cast<size_t>(by * nx + bx);
                samples[index]++;
                if (rc.occupied_set.count(encode_coord(x, y, layer))) counts[index]++;
            }
        }
        return json{
            {"bounds", json::array({xmin, xmax, ymin, ymax})},
            {"nx", nx},
            {"ny", ny},
            {"layer", layer},
            {"occupied_counts", counts},
            {"sample_counts", samples},
        };
    }

    json capture_snapshot(const BumpPair& pair, const RoutingCase& rc) const {
        const ExitTerminal& start = pair.first.exit_terminal;
        const ExitTerminal& end = pair.second.exit_terminal;
        return json{
            {"overview_density", capture_density(pair, rc)},
            {"start_local", capture_local_x2(2 * start.grid_coord_x,
                                               2 * start.grid_coord_y,
                                               start.layer, rc)},
            {"end_local", capture_local_x2(2 * end.grid_coord_x,
                                             2 * end.grid_coord_y,
                                             end.layer, rc)},
        };
    }

    json astar_neighbor_snapshot(const BumpPair& pair, const RoutingCase& rc) const {
        static constexpr int directions[10][3] = {
            {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0},
            {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
            {0, 0, 1}, {0, 0, -1},
        };
        const ExitTerminal& start = pair.first.exit_terminal;
        const ExitTerminal& target = pair.second.exit_terminal;
        json result = json::array();
        for (const auto& direction : directions) {
            const int x = start.grid_coord_x + direction[0];
            const int y = start.grid_coord_y + direction[1];
            const int layer = start.layer + direction[2];
            std::string reason = "free";
            if (x < 0 || x >= rc.grid_max_width || y < 0 || y >= rc.grid_max_height ||
                layer < 0 || layer >= rc.layer_count) {
                reason = "out_of_bounds";
            } else if (rc.occupied_set.count(encode_coord(x, y, layer)) &&
                       !(x == target.grid_coord_x && y == target.grid_coord_y &&
                         layer == target.layer)) {
                reason = "occupied_node";
            } else if (direction[0] != 0 && direction[1] != 0) {
                const float mid_x = start.grid_coord_x + 0.5f * direction[0];
                const float mid_y = start.grid_coord_y + 0.5f * direction[1];
                if (rc.occupied_set.count(encode_coord(mid_x, mid_y, layer))) {
                    reason = "occupied_diagonal_midpoint";
                }
            }
            result.push_back(json{
                {"dx", direction[0]}, {"dy", direction[1]}, {"dl", direction[2]},
                {"x", x}, {"y", y}, {"layer", layer}, {"state", reason},
            });
        }
        return result;
    }
};

}  // namespace

AStarNodeReroute::AStarNodeReroute(int x_, int y_, int layer_)
    : x(x_), y(y_), layer(layer_), g(0), h(0), f(0), parent(nullptr) {}

bool AStarNodeReroute::operator>(const AStarNodeReroute& other) const {
    return f > other.f;
}

// 删除 get_coord_node_key_reroute 和 get_node_key_reroute 函数

int grid_distance_reroute(int x1, int y1, int l1, int x2, int y2, int l2) {
    return std::max(abs(x1 - x2), abs(y1 - y2)) + abs(l1 - l2);
}

std::vector<AStarNodeReroute> get_neighbors_reroute(const AStarNodeReroute& current,
                                                    int max_width, int max_height, int layer_count,
                                                    const OccGrid& occupied_set,
                                                    int end_x, int end_y, int end_layer,
                                                    const std::unordered_set<uint64_t>* self_cells) {
    std::vector<AStarNodeReroute> neighbors;
    neighbors.reserve(10);
    static constexpr int dirs[10][3] = {
        {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0},
        {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
        {0, 0, 1}, {0, 0, -1}
    };

    for (const auto& dir : dirs) {
        int dx = dir[0];
        int dy = dir[1];
        int dl = dir[2];

        int new_x = current.x + dx;
        int new_y = current.y + dy;
        int new_layer = current.layer + dl;

        if (new_x >= 0 && new_x < max_width && new_y >= 0 && new_y < max_height && new_layer >= 0 && new_layer < layer_count) {
            uint64_t coord = encode_coord(new_x, new_y, new_layer);
            if (occupied_set.find(coord) != occupied_set.end()) {
                const bool own_cell = self_cells != nullptr &&
                                      self_cells->count(coord) != 0;
                if (!(new_x == end_x && new_y == end_y && new_layer == end_layer) &&
                    !own_cell) {
                    continue;
                }
            }
            if (dx != 0 && dy != 0) {
                if (current.parent != nullptr) {
                    if ((new_x == current.parent->x && current.y == current.parent->y) ||
                        (new_y == current.parent->y && current.x == current.parent->x)) {
                        continue;
                    }
                }
                float mid_x = static_cast<float>(current.x) + static_cast<float>(dx) * 0.5f;
                float mid_y = static_cast<float>(current.y) + static_cast<float>(dy) * 0.5f;
                uint64_t key = encode_coord(mid_x, mid_y, new_layer);
                const bool own_midpoint = self_cells != nullptr &&
                                          self_cells->count(key) != 0;
                if (occupied_set.find(key) != occupied_set.end() && !own_midpoint) {
                    continue;
                }
            } else {
                if (dx != 0) {
                    if (current.parent != nullptr) {
                        if (new_x == current.parent->x) {
                            continue;
                        }
                    }
                } else if (dy != 0) {
                    if (current.parent != nullptr) {
                        if (new_y == current.parent->y) {
                            continue;
                        }
                    }
                }
            }
            neighbors.emplace_back(new_x, new_y, new_layer);
        }
    }
    return neighbors;
}

std::vector<AStarNodeReroute> reconstruct_path_reroute(AStarNodeReroute* end_node,
                                                       OccGrid& occupied_set,
                                                       bool commit_path) {
    std::vector<AStarNodeReroute> path;
    AStarNodeReroute* current = end_node;
    AStarNodeReroute* last = end_node;
    while (current != nullptr) {
        path.push_back(*current);
        if (commit_path) {
            occupied_set.insert(
                encode_coord(current->x, current->y, current->layer));
        }
        if ((last->x != current->x) && (last->y != current->y)) {
            float mid_x = static_cast<float>(last->x + current->x) * 0.5f;
            float mid_y = static_cast<float>(last->y + current->y) * 0.5f;
            if (commit_path) {
                occupied_set.insert(
                    encode_coord(mid_x, mid_y, current->layer));
            }
        }
        last = current;
        current = current->parent;
    }
    return path;
}

bool reroute_path_conflicts(
    const std::vector<AStarNodeReroute>& path,
    const OccGrid& newly_committed) {
    for (size_t node_index = 0; node_index < path.size(); ++node_index) {
        const AStarNodeReroute& node = path[node_index];
        if (newly_committed.count(
                encode_coord(node.x, node.y, node.layer))) {
            return true;
        }
        if (node_index == 0) continue;
        const AStarNodeReroute& previous = path[node_index - 1];
        if (previous.x != node.x && previous.y != node.y &&
            newly_committed.count(encode_coord(
                static_cast<float>(previous.x + node.x) * 0.5f,
                static_cast<float>(previous.y + node.y) * 0.5f,
                node.layer))) {
            return true;
        }
    }
    return false;
}

void commit_reroute_path(
    const std::vector<AStarNodeReroute>& path,
    OccGrid& occupied_set,
    OccGrid& newly_committed) {
    for (size_t node_index = 0; node_index < path.size(); ++node_index) {
        const AStarNodeReroute& node = path[node_index];
        const uint64_t key = encode_coord(node.x, node.y, node.layer);
        occupied_set.insert(key);
        newly_committed.insert(key);
        if (node_index == 0) continue;
        const AStarNodeReroute& previous = path[node_index - 1];
        if (previous.x != node.x && previous.y != node.y) {
            const uint64_t midpoint_key = encode_coord(
                static_cast<float>(previous.x + node.x) * 0.5f,
                static_cast<float>(previous.y + node.y) * 0.5f,
                node.layer);
            occupied_set.insert(midpoint_key);
            newly_committed.insert(midpoint_key);
        }
    }
}

AStarResultReroute a_star_find_path_reroute(const BumpPair& bp,
                                            int max_width, int max_height, int layer_count,
                                            OccGrid& occupied_set,
                                            bool use_weighted_heuristic = false,
                                            long long node_budget = 0,
                                            const std::unordered_set<uint64_t>* self_cells = nullptr,
                                            bool commit_path = true) {
    const ExitTerminal& start_exit = bp.first.exit_terminal;
    const ExitTerminal& target_exit = bp.second.exit_terminal;
    int end_x = target_exit.grid_coord_x;
    int end_y = target_exit.grid_coord_y;
    int end_layer = target_exit.layer;

    std::priority_queue<AStarNodeReroute, std::vector<AStarNodeReroute>, std::greater<AStarNodeReroute>> open_list;
    static thread_local AStarStateGrid search_state;
    search_state.reset(max_width, max_height, layer_count);
    std::deque<AStarNodeReroute> parent_nodes;
    AStarResultReroute routing_result;

    AStarNodeReroute start_node(start_exit.grid_coord_x, start_exit.grid_coord_y, start_exit.layer);
    start_node.h = grid_distance_reroute(
        start_exit.grid_coord_x, start_exit.grid_coord_y, start_exit.layer,
        end_x, end_y, end_layer
    );
    start_node.f = use_weighted_heuristic
        ? 5 * start_node.g + 6 * start_node.h
        : start_node.g + start_node.h;
    open_list.push(start_node);
    search_state.set_best(
        start_node.x, start_node.y, start_node.layer, 0);

    while (!open_list.empty()) {
        AStarNodeReroute current = open_list.top();
        open_list.pop();
        routing_result.nodes_explored++;
        if (node_budget > 0 && routing_result.nodes_explored > node_budget) {
            routing_result.path = {};
            return routing_result;
        }

        if (current.x == end_x && current.y == end_y && current.layer == end_layer) {
            routing_result.path = reconstruct_path_reroute(
                &current, occupied_set, commit_path);
            return routing_result;
        }

        if (search_state.is_closed(current.x, current.y, current.layer) &&
            search_state.best_g(current.x, current.y, current.layer) < current.g) continue;
        search_state.close(current.x, current.y, current.layer);

        std::vector<AStarNodeReroute> neighbors = get_neighbors_reroute(current, max_width, max_height, layer_count,
                                                                        occupied_set, end_x, end_y, end_layer,
                                                                        self_cells);
        for (AStarNodeReroute& neighbor : neighbors) {
            if (search_state.is_closed(neighbor.x, neighbor.y, neighbor.layer)) continue;

            int tentative_g = current.g + 1;
            if (search_state.improve(neighbor.x, neighbor.y, neighbor.layer, tentative_g)) {
                parent_nodes.push_back(current);
                neighbor.parent = &parent_nodes.back();
                neighbor.g = tentative_g;
                neighbor.h = grid_distance_reroute(
                    neighbor.x, neighbor.y, neighbor.layer,
                    end_x, end_y, end_layer
                );
                neighbor.f = use_weighted_heuristic
                    ? 5 * neighbor.g + 6 * neighbor.h
                    : neighbor.g + neighbor.h;
                open_list.push(neighbor);
            }
        }
    }

    routing_result.path = {};
    return routing_result;
}

json coord_to_json_reroute(int x, int y, int layer) {
    std::string layer_name = (layer == 0) ? "Bottom" : "Top";
    return json::array({x, y, layer_name});
}

bool is_same_direction_reroute(int dx1, int dy1, int dl1,
                               int dx2, int dy2, int dl2) {
    return (dx1 == dx2) && (dy1 == dy2) && (dl1 == dl2);
}

std::vector<std::pair<ResultPoint, ResultPoint>> split_to_segments_reroute(const std::vector<AStarNodeReroute>& path) {
    std::vector<std::pair<ResultPoint, ResultPoint>> segments;
    if (path.size() < 2) {
        return segments;
    }

    ResultPoint current_start = {path[0].x, path[0].y, path[0].layer};
    int prev_dx = path[1].x - path[0].x;
    int prev_dy = path[1].y - path[0].y;
    int prev_dl = path[1].layer - path[0].layer;

    for (size_t i = 2; i < path.size(); ++i) {
        const AStarNodeReroute& prev_node = path[i-1];
        const AStarNodeReroute& curr_node = path[i];
        int curr_dx = curr_node.x - prev_node.x;
        int curr_dy = curr_node.y - prev_node.y;
        int curr_dl = curr_node.layer - prev_node.layer;

        if (!is_same_direction_reroute(prev_dx, prev_dy, prev_dl,
                                       curr_dx, curr_dy, curr_dl)) {
            ResultPoint new_node = {prev_node.x, prev_node.y, prev_node.layer};
            segments.emplace_back(current_start, new_node);
            current_start = new_node;
            prev_dx = curr_dx;
            prev_dy = curr_dy;
            prev_dl = curr_dl;
        }
    }

    const AStarNodeReroute& last_node = path.back();
    ResultPoint end_point = {last_node.x, last_node.y, last_node.layer};
    segments.emplace_back(current_start, end_point);
    return segments;
}

// 删除 get_node_key_2 函数，使用 encode_coord

ExitBumpPos getExitBumpBorder(const PairTeam &pt) {
    ExitBumpPos ebp;
    for (const BumpPair &bp : pt.member) {
        if (bp.first.exit_terminal.grid_coord_x == pt.start_left_bound_x) ebp.first.push_back(LEFT);
        else if (bp.first.exit_terminal.grid_coord_x == pt.start_right_bound_x) ebp.first.push_back(RIGHT);
        else if (bp.first.exit_terminal.grid_coord_y == pt.start_up_bound_y) ebp.first.push_back(DOWN);
        else if (bp.first.exit_terminal.grid_coord_y == pt.start_down_bound_y) ebp.first.push_back(UP);
        else IOUtils::terminateProgram("ERROR: Wrong Border Assignment!");

        if (bp.second.exit_terminal.grid_coord_x == pt.end_left_bound_x) ebp.second.push_back(LEFT);
        else if (bp.second.exit_terminal.grid_coord_x == pt.end_right_bound_x) ebp.second.push_back(RIGHT);
        else if (bp.second.exit_terminal.grid_coord_y == pt.end_up_bound_y) ebp.second.push_back(DOWN);
        else if (bp.second.exit_terminal.grid_coord_y == pt.end_down_bound_y) ebp.second.push_back(UP);
        else IOUtils::terminateProgram("ERROR: Wrong Border Assignment!");
    }
    return ebp;
}

bool canDirectRoute(const ResultPoint &a, const ResultPoint &b) {
    if (a.x == b.x || a.y == b.y || std::abs(a.y - b.y) == std::abs(a.x - b.x)) return true;
    return false;
}

bool canBendOneTime(ExitTerminal start, ExitPos start_pos, ExitTerminal end, ExitPos end_pos, bool &outStart, ResultPoint &cp) {
    if (!outStart) {
        std::swap(start, end);
        std::swap(start_pos, end_pos);
    }

    if (start_pos == end_pos) return false;
    if (start_pos == RIGHT && end_pos == LEFT && start.grid_coord_x > end.grid_coord_x) return false;
    if (start_pos == LEFT && end_pos == RIGHT && start.grid_coord_x < end.grid_coord_x) return false;
    if (start_pos == UP && end_pos == DOWN && start.grid_coord_y > end.grid_coord_y) return false;
    if (start_pos == DOWN && end_pos == UP && start.grid_coord_y < end.grid_coord_y) return false;

    int run_time = 0;
    while (++run_time < 3) {
        if (start_pos == UP) {
            cp.x = start.grid_coord_x;
            cp.y = end.grid_coord_y - std::abs(start.grid_coord_x - end.grid_coord_x);
            if (cp.y > start.grid_coord_y) return true;
        } else if (start_pos == RIGHT) {
            cp.y = start.grid_coord_y;
            cp.x = end.grid_coord_x - std::abs(start.grid_coord_y - end.grid_coord_y);
            if (cp.x > start.grid_coord_x) return true;
        } else if (start_pos == LEFT) {
            cp.y = start.grid_coord_y;
            cp.x = end.grid_coord_x + std::abs(start.grid_coord_y - end.grid_coord_y);
            if (cp.x < start.grid_coord_x) return true;
        } else {
            cp.x = start.grid_coord_x;
            cp.y = end.grid_coord_y + std::abs(start.grid_coord_x - end.grid_coord_x);
            if (cp.y < start.grid_coord_y) return true;
        }
        std::swap(start, end);
        std::swap(start_pos, end_pos);
        outStart = !outStart;
    }

    return false;
}

bool canBendTwoTime(const ExitTerminal &start, ExitPos start_pos, const ExitTerminal &end, ExitPos end_pos,
                    ResultPoint &cp1, ResultPoint &cp2, int bend_time, bool &no_space) {
    no_space = false;
    int abs_y_diff = std::abs(start.grid_coord_y - end.grid_coord_y);
    int abs_x_diff = std::abs(start.grid_coord_x - end.grid_coord_x);

    if (start_pos == LEFT && end_pos == RIGHT && abs_x_diff < abs_y_diff) {
        cp2.x = cp1.x = start.grid_coord_x - 2 * bend_time;
        if (cp1.x <= end.grid_coord_x) { no_space = true; return false; }
        int p = std::abs(start.grid_coord_y - end.grid_coord_y) / (start.grid_coord_y - end.grid_coord_y);
        cp1.y = start.grid_coord_y - p * 2 * bend_time;
        cp2.y = end.grid_coord_y + p * std::abs(end.grid_coord_x - cp1.x);
        return true;
    } else if (start_pos == RIGHT && end_pos == LEFT && abs_x_diff < abs_y_diff) {
        cp2.x = cp1.x = start.grid_coord_x + 2 * bend_time;
        if (cp1.x >= end.grid_coord_x) { no_space = true; return false; }
        int p = std::abs(start.grid_coord_y - end.grid_coord_y) / (start.grid_coord_y - end.grid_coord_y);
        cp1.y = start.grid_coord_y - p * 2 * bend_time;
        cp2.y = end.grid_coord_y + p * std::abs(end.grid_coord_x - cp1.x);
        return true;
    } else if (start_pos == UP && end_pos == DOWN && abs_x_diff > abs_y_diff) {
        int p = std::abs(start.grid_coord_x - end.grid_coord_x) / (start.grid_coord_x - end.grid_coord_x);
        cp1.y = cp2.y = start.grid_coord_y + 2 * bend_time;
        if (cp1.y >= end.grid_coord_y) { no_space = true; return false; }
        cp1.x = start.grid_coord_x - p * 2 * bend_time;
        cp2.x = end.grid_coord_x + p * std::abs(end.grid_coord_y - cp1.y);
        return true;
    } else if (start_pos == DOWN && end_pos == UP && abs_x_diff > abs_y_diff) {
        int p = std::abs(start.grid_coord_x - end.grid_coord_x) / (start.grid_coord_x - end.grid_coord_x);
        cp1.y = cp2.y = start.grid_coord_y - 2 * bend_time;
        if (cp1.y <= end.grid_coord_y) { no_space = true; return false; }
        cp1.x = start.grid_coord_x - p * 2 * bend_time;
        cp2.x = end.grid_coord_x + p * std::abs(end.grid_coord_y - cp1.y);
        return true;
    }

    return false;
}

bool update_occupied_grid(
    std::vector<Seg>& PerNetSeg,
    OccGrid& occupied_set,
    const OccGrid* base_occupied = nullptr) {

    auto is_occupied = [&](uint64_t key) {
        return occupied_set.count(key) ||
               (base_occupied != nullptr && base_occupied->count(key));
    };
    std::vector<uint64_t> inserting_keys;
    int seg_count = PerNetSeg.size();
    for (int seg_idx = 0; seg_idx < seg_count; ++seg_idx) {
        Seg &seg = PerNetSeg[seg_idx];
        ResultPoint &start = seg.first;
        ResultPoint &end = seg.second;
        if (start.layer != end.layer) {
            if (start.x != end.x || start.y != end.y) return false;
            const int first_layer = std::min(start.layer, end.layer);
            const int last_layer = std::max(start.layer, end.layer);
            for (int layer = first_layer; layer <= last_layer; ++layer) {
                const uint64_t key = encode_coord(start.x, start.y, layer);
                if (!is_occupied(key)) {
                    inserting_keys.push_back(key);
                } else if (!endpoint_overlap_allowed(
                               seg_idx, seg_count, key, start, end)) {
                    return false;
                }
            }
            continue;
        }
        int dx = end.x - start.x;
        int dy = end.y - start.y;
        int step_x = dx == 0 ? 0 : (dx > 0 ? 1 : -1);
        int step_y = dy == 0 ? 0 : (dy > 0 ? 1 : -1);
        int total_steps = std::max(abs(dx), abs(dy));
        int curr_x = start.x;
        int curr_y = start.y;
        if (dx != 0 && dy != 0) {
            for (int i = 0; i < total_steps; ++i) {
                uint64_t node_key = encode_coord(curr_x, curr_y, start.layer);
                if (!is_occupied(node_key)) {
                    inserting_keys.push_back(node_key);
                } else {
                    // 如果是当前线段的起点或终点，允许重复；否则冲突
                    if (!((seg_idx == 0 && node_key == encode_coord(start.x, start.y, start.layer)) ||
                          (seg_idx == seg_count - 1 && node_key == encode_coord(end.x, end.y, end.layer)))) {
                        return false;
                    }
                }
                float mid_x = static_cast<float>(curr_x) + static_cast<float>(step_x) * 0.5f;
                float mid_y = static_cast<float>(curr_y) + static_cast<float>(step_y) * 0.5f;
                uint64_t mid_key = encode_coord(mid_x, mid_y, start.layer);
                if (!is_occupied(mid_key)) {
                    inserting_keys.push_back(mid_key);
                } else {
                    if (!((seg_idx == 0 && mid_key == encode_coord(static_cast<float>(start.x), static_cast<float>(start.y), start.layer)) ||
                          (seg_idx == seg_count - 1 && mid_key == encode_coord(static_cast<float>(end.x), static_cast<float>(end.y), end.layer)))) {
                        return false;
                    }
                }
                curr_x += step_x;
                curr_y += step_y;
            }
            uint64_t last_node_key = encode_coord(end.x, end.y, start.layer);
            if (!is_occupied(last_node_key)) {
                inserting_keys.push_back(last_node_key);
            } else {
                if (!((seg_idx == 0 && last_node_key == encode_coord(start.x, start.y, start.layer)) ||
                      (seg_idx == seg_count - 1 && last_node_key == encode_coord(end.x, end.y, end.layer)))) {
                    return false;
                }
            }
        } else {
            for (int i = 0; i <= total_steps; ++i) {
                uint64_t node_key = encode_coord(curr_x, curr_y, start.layer);
                if (!is_occupied(node_key)) {
                    inserting_keys.push_back(node_key);
                } else {
                    if (!((seg_idx == 0 && node_key == encode_coord(start.x, start.y, start.layer)) ||
                          (seg_idx == seg_count - 1 && node_key == encode_coord(end.x, end.y, end.layer)))) {
                        return false;
                    }
                }
                curr_x += step_x;
                curr_y += step_y;
            }
        }
    }
    for (uint64_t key : inserting_keys) {
        occupied_set.insert(key);
    }
    return true;
}

bool update_occupied_set(std::vector<Seg> &PerNetSeg, RoutingCase& rc) {
    return update_occupied_grid(PerNetSeg, rc.occupied_set);
}

DirectRoutingResult direct_route(
    const std::vector<PairTeam>& vpt,
    RoutingCase& rc,
    const DRResultDict& escape_paths,
    int routing_threads) {
    const auto area_function_start = std::chrono::steady_clock::now();
    DirectRoutingResult res;
    DRResultDict ns;
    std::vector<PendingAreaNet> pending_nets;
    std::vector<std::string> final_failed_nets;
    size_t area_net_count = 0;
    AreaPatternTracer tracer;
    tracer.capture_area_start(vpt, rc);
    for (size_t team_index = 0; team_index < vpt.size(); ++team_index) {
        const PairTeam& pt = vpt[team_index];
        if (pt.member[0].first.exit_terminal == pt.member[0].second.exit_terminal) {
            for (const BumpPair & bp : pt.member) {
                std::vector<Seg> PerNetSeg;
                ns[bp.first.net_name] = PerNetSeg;
            }
            continue;
        }
        area_net_count += pt.member.size();

        ExitBumpPos ebp = getExitBumpBorder(pt);
        bool routeStart = true;
        int twobendtimes = 0;
        for (size_t i = 0; i < pt.member.size(); i++) {
            std::vector<Seg> PerNetSeg;
            const BumpPair bp = pt.member[i];
            const std::string& net_name = bp.first.net_name;
            bool need_reroute = false;
            std::string pattern_type = "none";
            ResultPoint rp1{bp.first.exit_terminal.grid_coord_x, bp.first.exit_terminal.grid_coord_y, bp.first.exit_terminal.layer};
            ResultPoint rp2{bp.second.exit_terminal.grid_coord_x, bp.second.exit_terminal.grid_coord_y, bp.second.exit_terminal.layer};
            if (canDirectRoute(rp1, rp2)) {
                pattern_type = "direct";
                PerNetSeg.emplace_back(rp1, rp2);
                ns[net_name] = PerNetSeg;
                twobendtimes = 0;
            } else {
                ResultPoint changePoint1{-1, -1, bp.first.exit_terminal.layer};
                ResultPoint changePoint2{-1, -1, bp.first.exit_terminal.layer};
                if (canBendOneTime(bp.first.exit_terminal, ebp.first[i], bp.second.exit_terminal, ebp.second[i], routeStart, changePoint1)) {
                    pattern_type = "one_bend";
                    PerNetSeg.emplace_back(rp1, changePoint1);
                    PerNetSeg.emplace_back(changePoint1, rp2);
                    ns[net_name] = PerNetSeg;
                    twobendtimes = 0;
                } else {
                    bool no_space = false;
                    if (canBendTwoTime(bp.first.exit_terminal, ebp.first[i], bp.second.exit_terminal, ebp.second[i], changePoint1, changePoint2, ++twobendtimes, no_space)) {
                        pattern_type = "two_bend";
                        PerNetSeg.emplace_back(rp1, changePoint1);
                        PerNetSeg.emplace_back(changePoint1, changePoint2);
                        PerNetSeg.emplace_back(changePoint2, rp2);
                    ns[net_name] = PerNetSeg;
                } else {
                    need_reroute = true;
                    ns[net_name] = {};
                }
            }

            if (need_reroute) {
                tracer.record_attempt(
                    bp, team_index, i, pattern_type, PerNetSeg,
                    "no_pattern_geometry", PatternConflict{}, rc);
                pending_nets.push_back({bp, {}, team_index, i});
                continue;
            }

            }  // close the } else { wrapper from canBendOneTime
            PatternConflict trace_conflict;
            if (tracer.wants(net_name)) {
                trace_conflict = find_first_pattern_conflict(PerNetSeg, rc);
                tracer.record_attempt(
                    bp, team_index, i, pattern_type, PerNetSeg,
                    trace_conflict.found ? "occupancy_conflict" : "committed",
                    trace_conflict, rc);
            }
            bool routing_success = update_occupied_set(PerNetSeg, rc);
            if (!routing_success) {
                pending_nets.push_back({bp, PerNetSeg, team_index, i});
            }
        }
    }
    const auto area_legacy_end = std::chrono::steady_clock::now();

    struct CandidateBundle {
        PendingAreaNet pending;
        std::vector<PatternCandidate> candidates;
    };
    std::vector<CandidateBundle> bundles(pending_nets.size());
    std::vector<size_t> generated_counts(pending_nets.size(), 0);
    constexpr size_t MAX_VIABLE_CANDIDATES_PER_LAYER = 32;
    const bool allow_layer_change =
        pending_nets.size() >= 64 &&
        pending_nets.size() * 10 >= std::max<size_t>(1, area_net_count);
    const int candidate_threads = std::max(
        1, std::min<int>(
            routing_threads <= 0 ? omp_get_max_threads() : routing_threads,
            static_cast<int>(std::max<size_t>(1, pending_nets.size()))));
    #pragma omp parallel for schedule(dynamic, 1) num_threads(candidate_threads)
    for (size_t pending_index = 0;
         pending_index < pending_nets.size(); ++pending_index) {
        const PendingAreaNet& pending = pending_nets[pending_index];
        std::vector<PatternCandidate> generated =
            generate_pattern_candidates(pending, rc, allow_layer_change);
        generated_counts[pending_index] = generated.size();
        CandidateBundle bundle{pending, {}};
        const size_t candidate_capacity = MAX_VIABLE_CANDIDATES_PER_LAYER *
            static_cast<size_t>(allow_layer_change ? rc.layer_count : 1);
        bundle.candidates.reserve(std::min(candidate_capacity, generated.size()));
        std::vector<size_t> candidates_per_layer(
            static_cast<size_t>(rc.layer_count), 0);
        for (PatternCandidate& candidate : generated) {
            int routing_layer = pending.pair.first.exit_terminal.layer;
            for (const Seg& segment : candidate.segments) {
                if (segment.first.layer == segment.second.layer &&
                    (segment.first.x != segment.second.x ||
                     segment.first.y != segment.second.y)) {
                    routing_layer = segment.first.layer;
                    break;
                }
            }
            size_t& layer_count =
                candidates_per_layer[static_cast<size_t>(routing_layer)];
            if (layer_count >= MAX_VIABLE_CANDIDATES_PER_LAYER) continue;
            if (find_first_pattern_conflict(candidate.segments, rc).found) continue;
            bundle.candidates.push_back(std::move(candidate));
            layer_count++;
        }
        bundles[pending_index] = std::move(bundle);
    }
    const size_t generated_candidates = std::accumulate(
        generated_counts.begin(), generated_counts.end(), size_t{0});
    const auto area_candidate_end = std::chrono::steady_clock::now();

    // Route the most constrained nets first. Successful legacy Patterns are
    // fixed, so this phase cannot perturb their geometry or quality.
    if (allow_layer_change) {
        std::stable_sort(bundles.begin(), bundles.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.candidates.size() != rhs.candidates.size()) {
            return lhs.candidates.size() < rhs.candidates.size();
        }
        const ResultPoint lhs_start{
            lhs.pending.pair.first.exit_terminal.grid_coord_x,
            lhs.pending.pair.first.exit_terminal.grid_coord_y,
            lhs.pending.pair.first.exit_terminal.layer};
        const ResultPoint lhs_end{
            lhs.pending.pair.second.exit_terminal.grid_coord_x,
            lhs.pending.pair.second.exit_terminal.grid_coord_y,
            lhs.pending.pair.second.exit_terminal.layer};
        const ResultPoint rhs_start{
            rhs.pending.pair.first.exit_terminal.grid_coord_x,
            rhs.pending.pair.first.exit_terminal.grid_coord_y,
            rhs.pending.pair.first.exit_terminal.layer};
        const ResultPoint rhs_end{
            rhs.pending.pair.second.exit_terminal.grid_coord_x,
            rhs.pending.pair.second.exit_terminal.grid_coord_y,
            rhs.pending.pair.second.exit_terminal.layer};
        const int lhs_span = std::max(
            std::abs(lhs_end.x - lhs_start.x),
            std::abs(lhs_end.y - lhs_start.y));
        const int rhs_span = std::max(
            std::abs(rhs_end.x - rhs_start.x),
            std::abs(rhs_end.y - rhs_start.y));
        if (lhs_span != rhs_span) return lhs_span > rhs_span;
        return lhs.pending.pair.first.net_name < rhs.pending.pair.first.net_name;
        });
    }

    std::vector<BumpPair> failed_nets;
    failed_nets.reserve(pending_nets.size());
    std::vector<uint8_t> pre_routed(bundles.size(), 0);
    size_t pre_astar_routed = 0;
    for (size_t bundle_index = 0; bundle_index < bundles.size(); ++bundle_index) {
        CandidateBundle& bundle = bundles[bundle_index];
        if (!allow_layer_change || !bundle.candidates.empty()) continue;
        std::unordered_set<uint64_t> self_cells = collect_net_escape_cells(
            bundle.pending.pair, escape_paths);
        AStarResultReroute pre_route = a_star_find_path_reroute(
            bundle.pending.pair,
            rc.grid_max_width, rc.grid_max_height, rc.layer_count,
            rc.occupied_set, true, 5000000LL, &self_cells);
        if (pre_route.path.empty()) continue;
        ns[bundle.pending.pair.first.net_name] =
            split_to_segments_reroute(pre_route.path);
        pre_routed[bundle_index] = 1;
        pre_astar_routed++;
    }
    const auto area_pre_astar_end = std::chrono::steady_clock::now();
    size_t zero_candidate_nets = 0;
    size_t viable_candidate_total = 0;
    size_t viable_candidate_max = 0;
    for (const CandidateBundle& bundle : bundles) {
        viable_candidate_total += bundle.candidates.size();
        viable_candidate_max = std::max(
            viable_candidate_max, bundle.candidates.size());
        if (bundle.candidates.empty()) zero_candidate_nets++;
    }

    const int allocation_trials = bundles.empty()
        ? 0 : (allow_layer_change ? 4 : 1);
    std::vector<int> best_selection(bundles.size(), -1);
    size_t best_routed = 0;
    int best_vias = std::numeric_limits<int>::max();
    double best_wirelength = std::numeric_limits<double>::infinity();
    std::vector<size_t> base_order(bundles.size());
    std::iota(base_order.begin(), base_order.end(), 0);
    auto name_hash = [](const std::string& name) {
        uint64_t hash = 1469598103934665603ULL;
        for (unsigned char byte : name) {
            hash ^= byte;
            hash *= 1099511628211ULL;
        }
        return hash;
    };
    std::vector<std::vector<int>> trial_selections(
        static_cast<size_t>(allocation_trials));
    std::vector<size_t> trial_routed(
        static_cast<size_t>(allocation_trials), 0);
    std::vector<int> trial_vias(
        static_cast<size_t>(allocation_trials), 0);
    std::vector<double> trial_wirelength(
        static_cast<size_t>(allocation_trials), 0.0);
    #pragma omp parallel for schedule(static) num_threads(candidate_threads)
    for (int trial = 0; trial < allocation_trials; ++trial) {
        OccGrid trial_occupancy;
        trial_occupancy.init(
            std::max(rc.grid_max_width, rc.grid_max_height) + 10000,
            rc.layer_count);
        std::vector<size_t> order = base_order;
        if (trial == 1) {
            std::reverse(order.begin(), order.end());
        } else if (trial >= 2) {
            std::stable_sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
                const uint64_t lhs_hash = name_hash(
                    bundles[lhs].pending.pair.first.net_name) ^
                    (0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(trial));
                const uint64_t rhs_hash = name_hash(
                    bundles[rhs].pending.pair.first.net_name) ^
                    (0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(trial));
                return lhs_hash < rhs_hash;
            });
        }
        std::vector<int> selection(bundles.size(), -1);
        size_t routed = 0;
        int vias = 0;
        double wirelength = 0.0;
        for (size_t bundle_index : order) {
            if (pre_routed[bundle_index]) continue;
            CandidateBundle& bundle = bundles[bundle_index];
            const size_t candidate_count = bundle.candidates.size();
            if (candidate_count == 0) continue;
            const size_t offset = trial == 0 ? 0 :
                static_cast<size_t>(name_hash(
                    bundle.pending.pair.first.net_name) + trial) % candidate_count;
            for (size_t attempt = 0; attempt < candidate_count; ++attempt) {
                const size_t candidate_index = (offset + attempt) % candidate_count;
                PatternCandidate& candidate = bundle.candidates[candidate_index];
                if (!update_occupied_grid(
                        candidate.segments, trial_occupancy,
                        &rc.occupied_set)) {
                    continue;
                }
                selection[bundle_index] = static_cast<int>(candidate_index);
                routed++;
                vias += candidate.vias;
                wirelength += candidate.wirelength;
                break;
            }
        }
        trial_selections[static_cast<size_t>(trial)] = std::move(selection);
        trial_routed[static_cast<size_t>(trial)] = routed;
        trial_vias[static_cast<size_t>(trial)] = vias;
        trial_wirelength[static_cast<size_t>(trial)] = wirelength;
    }
    for (int trial = 0; trial < allocation_trials; ++trial) {
        const size_t routed = trial_routed[static_cast<size_t>(trial)];
        const int vias = trial_vias[static_cast<size_t>(trial)];
        const double wirelength =
            trial_wirelength[static_cast<size_t>(trial)];
        const bool better = routed > best_routed ||
            (routed == best_routed && vias < best_vias) ||
            (routed == best_routed && vias == best_vias &&
             wirelength < best_wirelength);
        if (better) {
            best_routed = routed;
            best_vias = vias;
            best_wirelength = wirelength;
            best_selection = std::move(
                trial_selections[static_cast<size_t>(trial)]);
        }
    }
    const auto area_allocation_end = std::chrono::steady_clock::now();

    for (size_t bundle_index = 0; bundle_index < bundles.size(); ++bundle_index) {
        CandidateBundle& bundle = bundles[bundle_index];
        if (pre_routed[bundle_index]) continue;
        const int candidate_index = best_selection[bundle_index];
        if (candidate_index < 0) {
            failed_nets.push_back(bundle.pending.pair);
            continue;
        }
        PatternCandidate& candidate =
            bundle.candidates[static_cast<size_t>(candidate_index)];
        if (!update_occupied_set(candidate.segments, rc)) {
            failed_nets.push_back(bundle.pending.pair);
            continue;
        }
        ns[bundle.pending.pair.first.net_name] = candidate.segments;
    }
    const auto area_commit_end = std::chrono::steady_clock::now();
    const size_t multipattern_routed = bundles.size() - failed_nets.size();
    std::cout << "[AreaMultiPattern] pending " << pending_nets.size()
              << ", layer mode " << (allow_layer_change ? "adaptive" : "same-layer")
              << ", generated " << generated_candidates
              << ", viable " << viable_candidate_total
              << ", zero-candidate nets " << zero_candidate_nets
              << ", max/net " << viable_candidate_max
              << ", allocation trials " << allocation_trials
              << ", pre-routed zero-candidate " << pre_astar_routed
              << ", routed " << multipattern_routed
              << ", residual " << failed_nets.size() << std::endl;

    tracer.finish_pattern_stage(vpt, rc);
    size_t astar_probe_limit = 0;
    if (const char* probe_env = std::getenv("PIROUTER_AREA_ASTAR_PROBE_NETS")) {
        astar_probe_limit = static_cast<size_t>(std::strtoull(probe_env, nullptr, 10));
    }
    const size_t residual_count = astar_probe_limit > 0
        ? std::min(astar_probe_limit, failed_nets.size())
        : failed_nets.size();
    constexpr size_t MAX_TRANSACTIONAL_RESIDUAL_NETS = 32;
    if (astar_probe_limit > 0 ||
        failed_nets.size() > MAX_TRANSACTIONAL_RESIDUAL_NETS) {
        size_t astar_processed = 0;
        for (const BumpPair& bp : failed_nets) {
            std::unordered_set<uint64_t> self_cells;
            if (allow_layer_change) {
                self_cells = collect_net_escape_cells(bp, escape_paths);
            }
            const auto astar_start = std::chrono::steady_clock::now();
            AStarResultReroute reroute_result = a_star_find_path_reroute(
                bp, rc.grid_max_width, rc.grid_max_height,
                rc.layer_count, rc.occupied_set, allow_layer_change,
                allow_layer_change ? 5000000LL : 0LL,
                allow_layer_change ? &self_cells : nullptr,
                true);
            const double astar_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - astar_start).count();
            if (astar_seconds >= 0.05 || reroute_result.path.empty()) {
                std::cout << "[AreaAStar] net " << bp.first.net_name
                          << ", nodes " << reroute_result.nodes_explored
                          << ", elapsed " << astar_seconds << " s, "
                          << (reroute_result.path.empty()
                                  ? "failed" : "success")
                          << std::endl;
            }
            if (reroute_result.path.empty() &&
                !(bp.first.exit_terminal.grid_coord_x ==
                      bp.second.exit_terminal.grid_coord_x &&
                  bp.first.exit_terminal.grid_coord_y ==
                      bp.second.exit_terminal.grid_coord_y)) {
                final_failed_nets.push_back(bp.first.net_name);
                ns[bp.first.net_name] = {};
            } else {
                ns[bp.first.net_name] =
                    split_to_segments_reroute(reroute_result.path);
            }
            astar_processed++;
            if (astar_probe_limit > 0 &&
                astar_processed >= astar_probe_limit) {
                std::cout << "[AreaAStarProbe] stopped after "
                          << astar_processed << " residual nets"
                          << std::endl;
                std::exit(0);
            }
        }
        std::cout << "[AreaResidualParallel] nets " << residual_count
                  << ", threads 1, retries 0, mode sequential-dense"
                  << std::endl;
    } else {
    std::vector<AStarResultReroute> residual_results(residual_count);
    std::vector<double> residual_seconds(residual_count, 0.0);
    const int residual_threads = astar_probe_limit > 0
        ? 1
        : std::max(1, std::min<int>(
              candidate_threads,
              static_cast<int>(std::max<size_t>(1, residual_count))));
    size_t residual_retries = 0;
    const size_t RESIDUAL_TRANSACTION_BATCH_SIZE =
        std::max<size_t>(1, residual_count);
    for (size_t batch_begin = 0;
         batch_begin < residual_count;
         batch_begin += RESIDUAL_TRANSACTION_BATCH_SIZE) {
        const size_t batch_end = std::min(
            residual_count,
            batch_begin + RESIDUAL_TRANSACTION_BATCH_SIZE);
        #pragma omp parallel for schedule(dynamic, 1) num_threads(residual_threads)
        for (size_t failed_index = batch_begin;
             failed_index < batch_end; ++failed_index) {
        const BumpPair& bp = failed_nets[failed_index];
        std::unordered_set<uint64_t> self_cells;
        if (allow_layer_change) {
            self_cells = collect_net_escape_cells(bp, escape_paths);
        }
        const auto astar_start = std::chrono::steady_clock::now();
        residual_results[failed_index] = a_star_find_path_reroute(
            bp, rc.grid_max_width, rc.grid_max_height, rc.layer_count,
            rc.occupied_set, allow_layer_change,
            allow_layer_change ? 5000000LL : 0LL,
            allow_layer_change ? &self_cells : nullptr,
            false);
        residual_seconds[failed_index] = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - astar_start).count();
        }

        OccGrid committed_residual;
        committed_residual.init(
            std::max(rc.grid_max_width, rc.grid_max_height) + 10000,
            rc.layer_count);
        for (size_t failed_index = batch_begin;
             failed_index < batch_end; ++failed_index) {
        const BumpPair& bp = failed_nets[failed_index];
        AStarResultReroute reroute_result =
            std::move(residual_results[failed_index]);
        double astar_seconds = residual_seconds[failed_index];
        if (!reroute_result.path.empty() &&
            reroute_path_conflicts(
                reroute_result.path, committed_residual)) {
            residual_retries++;
            std::unordered_set<uint64_t> self_cells;
            if (allow_layer_change) {
                self_cells = collect_net_escape_cells(bp, escape_paths);
            }
            const auto retry_start = std::chrono::steady_clock::now();
            reroute_result = a_star_find_path_reroute(
                bp, rc.grid_max_width, rc.grid_max_height,
                rc.layer_count, rc.occupied_set, allow_layer_change,
                allow_layer_change ? 5000000LL : 0LL,
                allow_layer_change ? &self_cells : nullptr,
                true);
            astar_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - retry_start).count();
            for (const AStarNodeReroute& node : reroute_result.path) {
                committed_residual.insert(
                    encode_coord(node.x, node.y, node.layer));
            }
            for (size_t node_index = 1;
                 node_index < reroute_result.path.size(); ++node_index) {
                const AStarNodeReroute& previous =
                    reroute_result.path[node_index - 1];
                const AStarNodeReroute& node =
                    reroute_result.path[node_index];
                if (previous.x != node.x && previous.y != node.y) {
                    committed_residual.insert(encode_coord(
                        static_cast<float>(previous.x + node.x) * 0.5f,
                        static_cast<float>(previous.y + node.y) * 0.5f,
                        node.layer));
                }
            }
        } else {
            commit_reroute_path(
                reroute_result.path, rc.occupied_set,
                committed_residual);
        }
        if (astar_seconds >= 0.05 || reroute_result.path.empty()) {
            std::cout << "[AreaAStar] net " << bp.first.net_name
                      << ", nodes " << reroute_result.nodes_explored
                      << ", elapsed " << astar_seconds << " s, "
                      << (reroute_result.path.empty() ? "failed" : "success") << std::endl;
        }
        if (reroute_result.path.empty() && !(bp.first.exit_terminal.grid_coord_x == bp.second.exit_terminal.grid_coord_x && bp.first.exit_terminal.grid_coord_y == bp.second.exit_terminal.grid_coord_y)) {
            final_failed_nets.push_back(bp.first.net_name);
            ns[bp.first.net_name] = {};
        } else {
            std::vector<Seg> PerNetSeg = split_to_segments_reroute(reroute_result.path);
            ns[bp.first.net_name] = PerNetSeg;
        }
        if (astar_probe_limit > 0 &&
            failed_index + 1 >= astar_probe_limit) {
            std::cout << "[AreaAStarProbe] stopped after "
                      << failed_index + 1
                      << " residual nets" << std::endl;
            std::exit(0);
        }
        }
    }
    std::cout << "[AreaResidualParallel] nets " << residual_count
              << ", threads " << residual_threads
              << ", batch_size " << RESIDUAL_TRANSACTION_BATCH_SIZE
              << ", retries " << residual_retries << std::endl;
    }
    const auto area_residual_end = std::chrono::steady_clock::now();
    const auto phase_seconds = [](const auto& begin, const auto& end) {
        return std::chrono::duration<double>(end - begin).count();
    };
    std::cout << "[AreaPhaseTiming] legacy "
              << phase_seconds(area_function_start, area_legacy_end)
              << " s, candidate_generation "
              << phase_seconds(area_legacy_end, area_candidate_end)
              << " s, pre_astar "
              << phase_seconds(area_candidate_end, area_pre_astar_end)
              << " s, allocation "
              << phase_seconds(area_pre_astar_end, area_allocation_end)
              << " s, commit "
              << phase_seconds(area_allocation_end, area_commit_end)
              << " s, residual_astar "
              << phase_seconds(area_commit_end, area_residual_end)
              << " s, threads " << candidate_threads << std::endl;
    res.dr_result_dict = ns;
    res.final_failed_nets = final_failed_nets;
    return res;
}
