#include "layer_assignment.h"
#include "ioutils.h"
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <tuple>

namespace {

constexpr int TEAM_BOUNDARY_EXTEND = 5;
constexpr int MAX_ROUTING_LAYERS = 6;

uint64_t get_team_split_min_work() {
    // Ordinary large Teams use a conservative threshold. Intersecting Teams
    // have a separate threshold/quota because they cannot benefit from the
    // cheaper A/B half-Team split and otherwise remain serial long tails.
    constexpr uint64_t DEFAULT_TEAM_SPLIT_MIN_WORK = 250000;
    const char* raw_value = std::getenv("TEAM_SPLIT_MIN_WORK");
    if (raw_value == nullptr || raw_value[0] == '\0') {
        return DEFAULT_TEAM_SPLIT_MIN_WORK;
    }
    errno = 0;
    char* parse_end = nullptr;
    const unsigned long long parsed = std::strtoull(raw_value, &parse_end, 10);
    if (errno != 0 || raw_value[0] == '-' || parse_end == raw_value ||
        *parse_end != '\0') {
        std::cerr << "Invalid TEAM_SPLIT_MIN_WORK: " << raw_value << std::endl;
        std::exit(2);
    }
    return static_cast<uint64_t>(parsed);
}

uint64_t get_intersect_team_split_min_work() {
    constexpr uint64_t DEFAULT_INTERSECT_TEAM_SPLIT_MIN_WORK = 200000;
    const char* raw_value = std::getenv("TEAM_SPLIT_INTERSECT_MIN_WORK");
    if (raw_value == nullptr || raw_value[0] == '\0') {
        return DEFAULT_INTERSECT_TEAM_SPLIT_MIN_WORK;
    }
    errno = 0;
    char* parse_end = nullptr;
    const unsigned long long parsed = std::strtoull(raw_value, &parse_end, 10);
    if (errno != 0 || raw_value[0] == '-' || parse_end == raw_value ||
        *parse_end != '\0') {
        std::cerr << "Invalid TEAM_SPLIT_INTERSECT_MIN_WORK: " << raw_value
                  << std::endl;
        std::exit(2);
    }
    return static_cast<uint64_t>(parsed);
}

int get_max_accepted_team_splits() {
    constexpr int DEFAULT_MAX_ACCEPTED_TEAM_SPLITS = 3;
    const char* raw_value = std::getenv("TEAM_SPLIT_MAX_ACCEPTED");
    if (raw_value == nullptr || raw_value[0] == '\0') {
        return DEFAULT_MAX_ACCEPTED_TEAM_SPLITS;
    }
    errno = 0;
    char* parse_end = nullptr;
    const long parsed = std::strtol(raw_value, &parse_end, 10);
    if (errno != 0 || parse_end == raw_value || *parse_end != '\0' ||
        parsed < 0 || parsed > 16) {
        std::cerr << "Invalid TEAM_SPLIT_MAX_ACCEPTED: " << raw_value
                  << " (expected 0..16)" << std::endl;
        std::exit(2);
    }
    return static_cast<int>(parsed);
}

int get_max_accepted_intersect_team_splits() {
    constexpr int DEFAULT_MAX_ACCEPTED_INTERSECT_TEAM_SPLITS = 1;
    const char* raw_value =
        std::getenv("TEAM_SPLIT_MAX_INTERSECT_ACCEPTED");
    if (raw_value == nullptr || raw_value[0] == '\0') {
        return DEFAULT_MAX_ACCEPTED_INTERSECT_TEAM_SPLITS;
    }
    errno = 0;
    char* parse_end = nullptr;
    const long parsed = std::strtol(raw_value, &parse_end, 10);
    if (errno != 0 || parse_end == raw_value || *parse_end != '\0' ||
        parsed < 0 || parsed > 16) {
        std::cerr << "Invalid TEAM_SPLIT_MAX_INTERSECT_ACCEPTED: "
                  << raw_value << " (expected 0..16)" << std::endl;
        std::exit(2);
    }
    return static_cast<int>(parsed);
}

bool allow_disjoint_split_siblings_same_layer() {
    const char* raw_value = std::getenv("TEAM_SPLIT_SIBLING_LAYER_POLICY");
    if (raw_value == nullptr || raw_value[0] == '\0' ||
        std::string(raw_value) == "different") {
        return false;
    }
    if (std::string(raw_value) == "disjoint_same") return true;
    if (std::string(raw_value) == "preserve_parent") return false;
    std::cerr << "Invalid TEAM_SPLIT_SIBLING_LAYER_POLICY: " << raw_value
              << " (expected different, disjoint_same, or preserve_parent)"
              << std::endl;
    std::exit(2);
}

bool preserve_parent_split_layer() {
    const char* raw_value = std::getenv("TEAM_SPLIT_SIBLING_LAYER_POLICY");
    return raw_value != nullptr &&
           std::string(raw_value) == "preserve_parent";
}

bool endpoint_bounds_overlap(const PairTeam& team) {
    return !(team.start_left_bound_x > team.end_right_bound_x ||
             team.start_right_bound_x < team.end_left_bound_x ||
             team.start_up_bound_y > team.end_down_bound_y ||
             team.start_down_bound_y < team.end_up_bound_y);
}

void rebuild_team_geometry(
    PairTeam& team,
    int grid_width,
    int grid_height) {
    team.cnt = team.member.size();
    team.start_left_bound_x = team.end_left_bound_x =
        team.start_up_bound_y = team.end_up_bound_y = INF;
    team.start_right_bound_x = team.end_right_bound_x =
        team.start_down_bound_y = team.end_down_bound_y = -1;
    double sum_dx = 0.0;
    double sum_dy = 0.0;
    for (const BumpPair& pair : team.member) {
        const double dx = static_cast<double>(
            pair.second.grid_coord_x - pair.first.grid_coord_x);
        const double dy = static_cast<double>(
            pair.second.grid_coord_y - pair.first.grid_coord_y);
        const double pair_length = std::hypot(dx, dy);
        if (pair_length > 1e-9) {
            sum_dx += dx / pair_length;
            sum_dy += dy / pair_length;
        }
    }
    const double direction_length = std::hypot(sum_dx, sum_dy);
    if (direction_length > 1e-9) {
        team.aver_unit_vector = {
            sum_dx / direction_length, sum_dy / direction_length};
    }
    team.update_bounds();
    team.start_left_bound_x = std::max(
        0, team.start_left_bound_x - TEAM_BOUNDARY_EXTEND);
    team.start_right_bound_x = std::min(
        grid_width - 1, team.start_right_bound_x + TEAM_BOUNDARY_EXTEND);
    team.start_up_bound_y = std::max(
        0, team.start_up_bound_y - TEAM_BOUNDARY_EXTEND);
    team.start_down_bound_y = std::min(
        grid_height - 1, team.start_down_bound_y + TEAM_BOUNDARY_EXTEND);
    team.end_left_bound_x = std::max(
        0, team.end_left_bound_x - TEAM_BOUNDARY_EXTEND);
    team.end_right_bound_x = std::min(
        grid_width - 1, team.end_right_bound_x + TEAM_BOUNDARY_EXTEND);
    team.end_up_bound_y = std::max(
        0, team.end_up_bound_y - TEAM_BOUNDARY_EXTEND);
    team.end_down_bound_y = std::min(
        grid_height - 1, team.end_down_bound_y + TEAM_BOUNDARY_EXTEND);
}

uint64_t pair_work_proxy(const BumpPair& pair) {
    return static_cast<uint64_t>(std::max(
        std::abs(pair.second.grid_coord_x - pair.first.grid_coord_x),
        std::abs(pair.second.grid_coord_y - pair.first.grid_coord_y))) + 1;
}

uint64_t team_work_proxy(const PairTeam& team) {
    uint64_t work = 0;
    for (const BumpPair& pair : team.member) work += pair_work_proxy(pair);
    return work;
}

bool build_team_split_candidate(
    const PairTeam& parent,
    int grid_width,
    int grid_height,
    PairTeam& first_child,
    PairTeam& second_child) {
    const size_t member_count = parent.member.size();
    if (member_count < 4) return false;
    const size_t min_child_size = std::max<size_t>(
        2, std::min<size_t>(16, member_count / 4));
    if (member_count < 2 * min_child_size) return false;

    long long direction_x = 0;
    long long direction_y = 0;
    for (const BumpPair& pair : parent.member) {
        direction_x += pair.second.grid_coord_x - pair.first.grid_coord_x;
        direction_y += pair.second.grid_coord_y - pair.first.grid_coord_y;
    }

    auto order_key = [&](const BumpPair& pair, int mode) {
        const long long ax = pair.first.grid_coord_x;
        const long long ay = pair.first.grid_coord_y;
        const long long bx = pair.second.grid_coord_x;
        const long long by = pair.second.grid_coord_y;
        const long long midpoint_x2 = ax + bx;
        const long long midpoint_y2 = ay + by;
        switch (mode) {
            case 0: return ax;
            case 1: return ay;
            case 2: return bx;
            case 3: return by;
            case 4: return midpoint_x2;
            case 5: return midpoint_y2;
            case 6: return midpoint_x2 * direction_x + midpoint_y2 * direction_y;
            default: return -midpoint_x2 * direction_y + midpoint_y2 * direction_x;
        }
    };

    using CandidateScore = std::tuple<
        uint64_t, uint64_t, size_t, int, size_t>;
    CandidateScore best_score{
        std::numeric_limits<uint64_t>::max(),
        std::numeric_limits<uint64_t>::max(),
        std::numeric_limits<size_t>::max(), INT_MAX,
        std::numeric_limits<size_t>::max()};
    bool found = false;
    int best_mode = -1;
    size_t best_cut = 0;

    struct EndpointBounds {
        int start_left = INF;
        int start_right = -1;
        int start_top = INF;
        int start_bottom = -1;
        int end_left = INF;
        int end_right = -1;
        int end_top = INF;
        int end_bottom = -1;

        void add(const BumpPair& pair) {
            start_left = std::min(start_left, pair.first.grid_coord_x);
            start_right = std::max(start_right, pair.first.grid_coord_x);
            start_top = std::min(start_top, pair.first.grid_coord_y);
            start_bottom = std::max(start_bottom, pair.first.grid_coord_y);
            end_left = std::min(end_left, pair.second.grid_coord_x);
            end_right = std::max(end_right, pair.second.grid_coord_x);
            end_top = std::min(end_top, pair.second.grid_coord_y);
            end_bottom = std::max(end_bottom, pair.second.grid_coord_y);
        }

        void extend(int width, int height) {
            start_left = std::max(0, start_left - TEAM_BOUNDARY_EXTEND);
            start_right = std::min(
                width - 1, start_right + TEAM_BOUNDARY_EXTEND);
            start_top = std::max(0, start_top - TEAM_BOUNDARY_EXTEND);
            start_bottom = std::min(
                height - 1, start_bottom + TEAM_BOUNDARY_EXTEND);
            end_left = std::max(0, end_left - TEAM_BOUNDARY_EXTEND);
            end_right = std::min(
                width - 1, end_right + TEAM_BOUNDARY_EXTEND);
            end_top = std::max(0, end_top - TEAM_BOUNDARY_EXTEND);
            end_bottom = std::min(
                height - 1, end_bottom + TEAM_BOUNDARY_EXTEND);
        }
    };

    auto endpoint_overlap = [](const EndpointBounds& bounds) {
        return !(bounds.start_left > bounds.end_right ||
                 bounds.start_right < bounds.end_left ||
                 bounds.start_top > bounds.end_bottom ||
                 bounds.start_bottom < bounds.end_top);
    };
    auto sibling_overlap = [](
        const EndpointBounds& first,
        const EndpointBounds& second,
        bool start_side) {
        const int first_left = start_side
            ? first.start_left : first.end_left;
        const int first_right = start_side
            ? first.start_right : first.end_right;
        const int first_top = start_side
            ? first.start_top : first.end_top;
        const int first_bottom = start_side
            ? first.start_bottom : first.end_bottom;
        const int second_left = start_side
            ? second.start_left : second.end_left;
        const int second_right = start_side
            ? second.start_right : second.end_right;
        const int second_top = start_side
            ? second.start_top : second.end_top;
        const int second_bottom = start_side
            ? second.start_bottom : second.end_bottom;
        return !(first_left > second_right ||
                 first_right < second_left ||
                 first_top > second_bottom ||
                 first_bottom < second_top);
    };

    for (int mode = 0; mode < 8; ++mode) {
        std::vector<size_t> order(member_count);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
            const long long lhs_key = order_key(parent.member[lhs], mode);
            const long long rhs_key = order_key(parent.member[rhs], mode);
            if (lhs_key != rhs_key) return lhs_key < rhs_key;
            return parent.member[lhs].first.net_name <
                   parent.member[rhs].first.net_name;
        });

        std::vector<uint64_t> prefix_work(member_count + 1, 0);
        for (size_t i = 0; i < member_count; ++i) {
            prefix_work[i + 1] = prefix_work[i] +
                pair_work_proxy(parent.member[order[i]]);
        }
        std::vector<EndpointBounds> prefix_bounds(member_count + 1);
        std::vector<EndpointBounds> suffix_bounds(member_count + 1);
        for (size_t i = 0; i < member_count; ++i) {
            prefix_bounds[i + 1] = prefix_bounds[i];
            prefix_bounds[i + 1].add(parent.member[order[i]]);
            const size_t reverse_index = member_count - 1 - i;
            suffix_bounds[reverse_index] = suffix_bounds[reverse_index + 1];
            suffix_bounds[reverse_index].add(
                parent.member[order[reverse_index]]);
        }
        const uint64_t total_work = prefix_work.back();
        for (size_t cut = min_child_size;
             cut + min_child_size <= member_count; ++cut) {
            EndpointBounds candidate_a_bounds = prefix_bounds[cut];
            EndpointBounds candidate_b_bounds = suffix_bounds[cut];
            candidate_a_bounds.extend(grid_width, grid_height);
            candidate_b_bounds.extend(grid_width, grid_height);
            if (endpoint_overlap(candidate_a_bounds) ||
                endpoint_overlap(candidate_b_bounds)) {
                continue;
            }
            // ExitBuilder reserves exit coordinates in a 2-D set shared by all
            // layers. Keep sibling A boxes and sibling B boxes disjoint so the
            // split does not manufacture exit collisions before Escape starts.
            if (sibling_overlap(
                    candidate_a_bounds, candidate_b_bounds, true) ||
                sibling_overlap(
                    candidate_a_bounds, candidate_b_bounds, false)) {
                continue;
            }

            const uint64_t first_work = prefix_work[cut];
            const uint64_t second_work = total_work - first_work;
            const CandidateScore score{
                std::max(first_work, second_work),
                first_work > second_work
                    ? first_work - second_work : second_work - first_work,
                cut > member_count - cut
                    ? cut - (member_count - cut) : (member_count - cut) - cut,
                mode, cut};
            if (!found || score < best_score) {
                found = true;
                best_score = score;
                best_mode = mode;
                best_cut = cut;
            }
        }
    }
    if (found) {
        std::vector<size_t> order(member_count);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
            const long long lhs_key = order_key(parent.member[lhs], best_mode);
            const long long rhs_key = order_key(parent.member[rhs], best_mode);
            if (lhs_key != rhs_key) return lhs_key < rhs_key;
            return parent.member[lhs].first.net_name <
                   parent.member[rhs].first.net_name;
        });
        first_child = parent;
        second_child = parent;
        first_child.member.clear();
        second_child.member.clear();
        first_child.member.reserve(best_cut);
        second_child.member.reserve(member_count - best_cut);
        for (size_t i = 0; i < member_count; ++i) {
            (i < best_cut ? first_child.member : second_child.member)
                .push_back(parent.member[order[i]]);
        }
        rebuild_team_geometry(first_child, grid_width, grid_height);
        rebuild_team_geometry(second_child, grid_width, grid_height);
    }
    return found;
}

int count_same_layer_conflicts(
    const std::vector<std::vector<int>>& graph,
    const std::vector<int>& assignment) {
    if (assignment.size() != graph.size()) return INT_MAX;
    int conflicts = 0;
    for (int i = 0; i < static_cast<int>(graph.size()); ++i) {
        if (assignment[i] < 0 || assignment[i] >= MAX_ROUTING_LAYERS) {
            return INT_MAX;
        }
        for (int neighbor : graph[i]) {
            if (i < neighbor && assignment[i] == assignment[neighbor]) {
                conflicts++;
            }
        }
    }
    return conflicts;
}

int used_layer_count(const std::vector<int>& assignment) {
    int max_layer = -1;
    for (int layer : assignment) max_layer = std::max(max_layer, layer);
    return max_layer + 1;
}

}  // namespace

// Calculate the occupancy area of a PairTeam
BusOccupancy calculate_occupancy(
    int start_left_bound_x,
    int start_right_bound_x,
    int start_up_bound_y,
    int start_down_bound_y,
    int end_left_bound_x,
    int end_right_bound_x,
    int end_up_bound_y,
    int end_down_bound_y
) {
    BusOccupancy occ;

    // 注释部分不改，但编号从0开始（若需要可调整，此处省略）
    occ.x_min = std::min(start_left_bound_x, end_left_bound_x);
    occ.x_max = std::max(start_right_bound_x, end_right_bound_x);
    occ.y_min = std::min(start_up_bound_y, end_up_bound_y);
    occ.y_max = std::max(start_down_bound_y, end_down_bound_y);

    std::vector<intPoint> start_rect_points = {
        {start_left_bound_x, start_up_bound_y},
        {start_right_bound_x, start_up_bound_y},
        {start_left_bound_x, start_down_bound_y},
        {start_right_bound_x, start_down_bound_y}
    };

    std::vector<intPoint> end_rect_points = {
        {end_left_bound_x, end_up_bound_y},
        {end_right_bound_x, end_up_bound_y},
        {end_left_bound_x, end_down_bound_y},
        {end_right_bound_x, end_down_bound_y}
    };

    std::unordered_map<int, int> start_corner_cuts;

    for (const intPoint& p : start_rect_points) {
        if (p.second != occ.y_min && p.second != occ.y_max) {
            if (p.first == occ.x_min) {
                if (p.second == start_up_bound_y) {
                    start_corner_cuts[1] = p.second - occ.y_min;
                } else { 
                    start_corner_cuts[2] = occ.y_max - p.second;
                }
            } else if (p.first == occ.x_max) {
                if (p.second == start_up_bound_y) {
                    start_corner_cuts[4] = p.second - occ.y_min;
                } else {
                    start_corner_cuts[3] = occ.y_max - p.second;
                }
            }
        } else if (p.first != occ.x_min && p.first != occ.x_max) {
            if (p.second == occ.y_min) {
                if (p.first == start_left_bound_x) {
                    start_corner_cuts[1] = p.first - occ.x_min;
                } else {
                    start_corner_cuts[4] = occ.x_max - p.first;
                }
            } else if (p.second == occ.y_max) {
                if (p.first == start_left_bound_x) {
                    start_corner_cuts[2] = p.first - occ.x_min;
                } else {
                    start_corner_cuts[3] = occ.x_max - p.first;
                }
            }
        }
    }

    for (const intPoint& p : end_rect_points) {
        if (p.second != occ.y_min && p.second != occ.y_max) {
            if (p.first == occ.x_min) {
                int update_position;
                int update_size;
                if (p.second == end_up_bound_y) {
                    update_position = 1;
                    update_size = p.second - occ.y_min;
                } else {
                    update_position = 2;
                    update_size = occ.y_max - p.second;
                }
                if (start_corner_cuts.find(update_position) != start_corner_cuts.end()) {
                    occ.corner_cuts[update_position] = std::min(start_corner_cuts[update_position], update_size);
                }
            } else if (p.first == occ.x_max) {
                int update_position;
                int update_size;
                if (p.second == end_up_bound_y) {
                    update_position = 4;
                    update_size = p.second - occ.y_min;
                } else {
                    update_position = 3;
                    update_size = occ.y_max - p.second;
                }
                if (start_corner_cuts.find(update_position) != start_corner_cuts.end()) {
                    occ.corner_cuts[update_position] = std::min(start_corner_cuts[update_position], update_size);
                }
            }
        } else if (p.first != occ.x_min && p.first != occ.x_max) {
            if (p.second == occ.y_min) {
                int update_position;
                int update_size;
                if (p.first == end_left_bound_x) {
                    update_position = 1;
                    update_size = p.first - occ.x_min;
                } else {
                    update_position = 4;
                    update_size = occ.x_max - p.first;
                }
                if (start_corner_cuts.find(update_position) != start_corner_cuts.end()) {
                    occ.corner_cuts[update_position] = std::min(start_corner_cuts[update_position], update_size);
                }
            } else if (p.second == occ.y_max) {
                int update_position;
                int update_size;
                if (p.first == end_left_bound_x) {
                    update_position = 2;
                    update_size = p.first - occ.x_min;
                } else {
                    update_position = 3;
                    update_size = occ.x_max - p.first;
                }
                if (start_corner_cuts.find(update_position) != start_corner_cuts.end()) {
                    occ.corner_cuts[update_position] = std::min(start_corner_cuts[update_position], update_size);
                }
            }
        }
    }

    return occ;
}

// Generate vertices of the convex polygon for the rectangle with cut corners (in clockwise order)
std::vector<Point> generate_cut_rect_vertices(const BusOccupancy& occ) {
    std::vector<Point> vertices;
    double x1 = static_cast<double>(occ.x_min);
    double x2 = static_cast<double>(occ.x_max);
    double y1 = static_cast<double>(occ.y_min);
    double y2 = static_cast<double>(occ.y_max);

    for (int corner_pos = 1;  corner_pos <= 4; corner_pos++) {
        if (occ.corner_cuts.count(corner_pos) == 0) {
            switch (corner_pos) {
                case 1: vertices.emplace_back(x1, y1); break;
                case 2: vertices.emplace_back(x1, y2); break;
                case 3: vertices.emplace_back(x2, y2); break;
                case 4: vertices.emplace_back(x2, y1); break;
            }
        } else {
            double s = static_cast<double>(occ.corner_cuts.at(corner_pos));
            switch (corner_pos) {
                case 1:
                    vertices.emplace_back(x1 + s, y1);
                    vertices.emplace_back(x1, y1 + s);
                    break;
                case 2:
                    vertices.emplace_back(x1, y2 - s);
                    vertices.emplace_back(x1 + s, y2);
                    break;
                case 3:
                    vertices.emplace_back(x2 - s, y2);
                    vertices.emplace_back(x2, y2 - s);
                    break;
                case 4:
                    vertices.emplace_back(x2, y1 + s);
                    vertices.emplace_back(x2 - s, y1);
                    break;
            }
        }
    }
    return vertices;
}

// Projection of a convex polygon on a specified axis
std::pair<double, double> LayerAssigner::project_convex_polygon(const std::vector<Point>& poly, const Point& axis) const {
    double min_proj = dot(poly[0], axis);
    double max_proj = min_proj;
    for (size_t i = 1; i < poly.size(); i++) {
        double proj = dot(poly[i], axis);
        min_proj = std::min(min_proj, proj);
        max_proj = std::max(max_proj, proj);
    }
    return {min_proj, max_proj};
}

// Separating Axis Theorem for determining if two convex polygons intersect
bool LayerAssigner::convex_polygons_intersect(const std::vector<Point>& polyA, const std::vector<Point>& polyB) const {
    // Extract all edge normals of the two polygons (candidates for separating axes)
    auto getEdgeNormals = [&](const std::vector<Point>& poly) -> std::vector<Point> {
        std::vector<Point> normals;
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            Point edge = poly[(i+1)%n] - poly[i];
            Point normal = {-edge.second, edge.first}; // rotate 90 degrees
            normal = normalize(normal);
            normals.push_back(normal);
        }
        return normals;
    };

    std::vector<Point> normalsA = getEdgeNormals(polyA);
    std::vector<Point> normalsB = getEdgeNormals(polyB);

    // Check all separating axes: if there exists a non-overlapping projection, the polygons do not intersect
    auto checkAxis = [&](const Point& axis) -> bool {
        std::pair<double, double> projA = project_convex_polygon(polyA, axis);
        std::pair<double, double> projB = project_convex_polygon(polyB, axis);
        return !(projA.second < projB.first - 1e-6 || projB.second < projA.first - 1e-6);
    };

    for (const auto& axis : normalsA) {
        if (!checkAxis(axis)) return false;
    }
    for (const auto& axis : normalsB) {
        if (!checkAxis(axis)) return false;
    }

    return true;
}

// Bus conflict detection
bool LayerAssigner::bus_conflict(const std::vector<Point>& polyA, const std::vector<Point>& polyB) const {
    return convex_polygons_intersect(polyA, polyB);
}

// Build conflict graph based on bus conflicts
std::vector<std::vector<int>> LayerAssigner::build_conflict_graph(const std::vector<PairTeam>& vpt) const {
    int n = vpt.size();
    std::vector<std::vector<int>> conflict_graph(n);
    std::vector<std::vector<Point>> all_resources;

    int count = 0;
    for (const auto& pt : vpt) {
        // 输出时编号改为从0开始，注释中保持一致性
        // std::cout << "Processing PairTeam " << count << "/" << n << std::endl;
        BusOccupancy occ = calculate_occupancy(
            pt.start_left_bound_x,
            pt.start_right_bound_x,
            pt.start_up_bound_y,
            pt.start_down_bound_y,
            pt.end_left_bound_x,
            pt.end_right_bound_x,
            pt.end_up_bound_y,
            pt.end_down_bound_y
        );
        // 注释部分不变
        std::vector<Point> r_resource = generate_cut_rect_vertices(occ);
        all_resources.push_back(r_resource);
        count++;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const bool split_siblings =
                vpt[i].split_family_id >= 0 &&
                vpt[i].split_family_id == vpt[j].split_family_id;
            const bool geometric_conflict =
                bus_conflict(all_resources[i], all_resources[j]);
            const bool force_sibling_conflict =
                split_siblings &&
                !allow_disjoint_split_siblings_same_layer();
            if (force_sibling_conflict || geometric_conflict) {
                conflict_graph[i].push_back(j);
                conflict_graph[j].push_back(i);
                // 输出时编号从0开始
                // std::cout << "PairTeam " << i << " and PairTeam " << j << " Conflict" << std::endl;
            }
        }
    }

    return conflict_graph;
}

// Backtracking for optimal layer assignment
void LayerAssigner::backtrack(
    int current_idx,
    const std::vector<std::vector<int>>& conflict_graph,
    std::vector<int>& current_assignment,
    int current_max_layer,
    std::vector<int>& best_assignment,
    int& best_total_layers
) const {
    int n = conflict_graph.size();
    // all teams assigned, check if better than best
    if (current_idx == n) {
        int current_total = current_max_layer + 1;
        if (current_total < best_total_layers) {
            best_total_layers = current_total;
            best_assignment = current_assignment;
        }
        return;
    }
    // Pruning: current layer has reached or exceeded best layer, no need to continue
    if (current_max_layer + 1 >= best_total_layers) return;
    // Collect used layers from conflict nodes
    std::unordered_set<int> used_layers;
    for (int conflict_idx : conflict_graph[current_idx]) {
        if (current_assignment[conflict_idx] != -1) {
            used_layers.insert(current_assignment[conflict_idx]);
        }
    }
    // Try all possible layer numbers
    for (int layer = 0; layer < best_total_layers; ++layer) {
        if (used_layers.count(layer)) continue; // Skip used layers
        current_assignment[current_idx] = layer;
        backtrack(current_idx + 1, conflict_graph, current_assignment, 
                  std::max(current_max_layer, layer), best_assignment, best_total_layers);
        current_assignment[current_idx] = -1; // Backtrack: undo assignment
    }
}

// Get optimal layer assignment (backtracking)
LayerAssignmentResult LayerAssigner::assign_optimal_layers(const std::vector<PairTeam>& vpt) {
    LayerAssignmentResult result;
    int n = vpt.size();
    if (n == 0) {
        result.total_layers = 0;
        result.is_optimal = true;
        return result;
    }
    // 日志中编号改为从0开始
    // IOUtils::print("Starting optimal layer assignment (backtracking)...\n");
    // IOUtils::print("Number of PairTeams: " + std::to_string(n) + "\n");
    // Build conflict graph
    auto conflict_graph = build_conflict_graph(vpt);
    // Initialize backtracking parameters
    std::vector<int> current_assignment(n, -1);  // -1 means unassigned
    std::vector<int> best_assignment(n, -1);
    int best_total_layers = INT_MAX;
    // Start backtracking
    backtrack(0, conflict_graph, current_assignment, -1, best_assignment, best_total_layers);
    result.team_layer = best_assignment;
    result.total_layers = best_total_layers;
    result.is_optimal = true;
    // IOUtils::print("Optimal layer assignment completed. Total layers needed: " + std::to_string(best_total_layers) + "\n");
    return result;
}

// Get approximate optimal solution (greedy)
LayerAssignmentResult LayerAssigner::assign_approximate_layers(const std::vector<PairTeam>& vpt) {
    LayerAssignmentResult result;
    int n = vpt.size();
    if (n == 0) {
        result.total_layers = 0;
        result.is_optimal = false;
        return result;
    }
    // 日志中编号改为从0开始
    // IOUtils::print("Starting approximate layer assignment (greedy)...\n");
    // IOUtils::print("Number of PairTeams: " + std::to_string(n) + "\n");
    // Build conflict graph
    auto conflict_graph = build_conflict_graph(vpt);
    // Sort by conflict degree (prioritize teams with more conflicts)
    std::vector<int> team_indices(n);
    std::iota(team_indices.begin(), team_indices.end(), 0);
    std::sort(team_indices.begin(), team_indices.end(), [&](int a, int b) {
        return conflict_graph[a].size() > conflict_graph[b].size();
    });
    // Greedily assign layer numbers
    // 策略：优先选"最小可用层号"（保持原有紧凑分配，不打乱本来就分得好的 case）；
    // 仅当该层已超载（net 数 > CAP）时，才跳到下一个未超载且非冲突的层，分流拥挤层。
    result.team_layer.resize(n, -1);
    int max_layer = 0;
    const long long PER_LAYER_CAP = 2500;  // 每层 net 软上限，仅超载才分流
    std::vector<long long> layer_load;
    for (int idx : team_indices) {
        std::unordered_set<int> used_layers;
        for (int conflict_idx : conflict_graph[idx]) {
            if (result.team_layer[conflict_idx] != -1) {
                used_layers.insert(result.team_layer[conflict_idx]);
            }
        }
        int team_size = (int)vpt[idx].member.size();
        // 从最小层号开始找：非冲突 且 (未超载 或 已到最后)
        int assigned_layer = 0;
        while (true) {
            if (used_layers.count(assigned_layer)) { assigned_layer++; continue; }
            long long load = (assigned_layer < (int)layer_load.size()) ? layer_load[assigned_layer] : 0;
            if (load > PER_LAYER_CAP) { assigned_layer++; continue; }  // 超载 → 换层分流
            break;
        }
        result.team_layer[idx] = assigned_layer;
        if (assigned_layer >= (int)layer_load.size()) layer_load.resize(assigned_layer + 1, 0);
        layer_load[assigned_layer] += team_size;
        max_layer = std::max(max_layer, assigned_layer);
    }
    result.total_layers = max_layer + 1;
    result.is_optimal = false;
    // IOUtils::print("Approximate layer assignment completed. Total layers needed: " + std::to_string(result.total_layers) + "\n");
    return result;
}

LayerAssignmentResult LayerAssigner::assign_bounded_layers(
    const std::vector<PairTeam>& vpt,
    const std::vector<std::vector<int>>& conflict_graph,
    const LayerAssignmentResult& baseline,
    int max_layers
) const {
    LayerAssignmentResult result;
    const int n = static_cast<int>(vpt.size());
    result.team_layer.assign(n, -1);
    result.total_layers = max_layers;
    result.is_optimal = false;
    if (n == 0 || max_layers <= 0) {
        result.total_layers = 0;
        return result;
    }

    std::vector<long long> team_size(n, 1);
    for (int i = 0; i < n; ++i) {
        team_size[i] = std::max<long long>(1, static_cast<long long>(vpt[i].member.size()));
    }

    struct StrictSolver {
        const std::vector<std::vector<int>>& graph;
        const std::vector<long long>& sizes;
        const std::vector<int>& baseline_layer;
        int layer_count;
        long long node_limit;
        std::chrono::steady_clock::time_point deadline;
        std::vector<int> color;
        std::vector<int> best_color;
        std::vector<long long> loads;
        long long nodes = 0;
        bool aborted = false;

        int select_team() const {
            int selected = -1;
            int best_saturation = -1;
            int best_degree = -1;
            long long best_size = -1;
            for (int i = 0; i < static_cast<int>(color.size()); ++i) {
                if (color[i] >= 0) continue;
                bool used[16] = {};
                int saturation = 0;
                for (int neighbor : graph[i]) {
                    int c = color[neighbor];
                    if (c >= 0 && !used[c]) {
                        used[c] = true;
                        saturation++;
                    }
                }
                int degree = static_cast<int>(graph[i].size());
                if (saturation > best_saturation ||
                    (saturation == best_saturation && degree > best_degree) ||
                    (saturation == best_saturation && degree == best_degree && sizes[i] > best_size) ||
                    (saturation == best_saturation && degree == best_degree && sizes[i] == best_size &&
                     (selected < 0 || i < selected))) {
                    selected = i;
                    best_saturation = saturation;
                    best_degree = degree;
                    best_size = sizes[i];
                }
            }
            return selected;
        }

        bool search() {
            if (++nodes > node_limit || std::chrono::steady_clock::now() >= deadline) {
                aborted = true;
                return false;
            }
            int idx = select_team();
            if (idx < 0) {
                best_color = color;
                return true;
            }

            bool forbidden[16] = {};
            for (int neighbor : graph[idx]) {
                int c = color[neighbor];
                if (c >= 0) forbidden[c] = true;
            }
            std::vector<int> candidates;
            for (int layer = 0; layer < layer_count; ++layer) {
                if (!forbidden[layer]) candidates.push_back(layer);
            }
            int preferred = baseline_layer[idx] < layer_count ? baseline_layer[idx] : layer_count - 1;
            std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
                bool a_preferred = a == preferred;
                bool b_preferred = b == preferred;
                if (a_preferred != b_preferred) return a_preferred;
                if (loads[a] != loads[b]) return loads[a] < loads[b];
                int da = std::abs(a - preferred);
                int db = std::abs(b - preferred);
                if (da != db) return da < db;
                return a < b;
            });
            for (int layer : candidates) {
                color[idx] = layer;
                loads[layer] += sizes[idx];
                if (search()) return true;
                loads[layer] -= sizes[idx];
                color[idx] = -1;
                if (aborted) return false;
            }
            return false;
        }
    };

    const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(250);
    const long long strict_node_limit = 2000000;
    std::vector<int> baseline_layer = baseline.team_layer;

    // First keep all teams already on legal layers fixed and insert only overflow teams.
    StrictSolver insertion_solver{
        conflict_graph, team_size, baseline_layer, max_layers,
        strict_node_limit / 4, deadline,
        std::vector<int>(n, -1), {}, std::vector<long long>(max_layers, 0)
    };
    for (int i = 0; i < n; ++i) {
        if (baseline_layer[i] >= 0 && baseline_layer[i] < max_layers) {
            insertion_solver.color[i] = baseline_layer[i];
            insertion_solver.loads[baseline_layer[i]] += team_size[i];
        }
    }
    bool strict_found = insertion_solver.search();
    long long strict_nodes = insertion_solver.nodes;
    bool preserved_existing = strict_found;

    // If overflow-only insertion is impossible, allow DSATUR to recolor all teams.
    if (!strict_found && std::chrono::steady_clock::now() < deadline) {
        StrictSolver global_solver{
            conflict_graph, team_size, baseline_layer, max_layers,
            strict_node_limit - std::min(strict_node_limit, strict_nodes), deadline,
            std::vector<int>(n, -1), {}, std::vector<long long>(max_layers, 0)
        };
        strict_found = global_solver.search();
        strict_nodes += global_solver.nodes;
        if (strict_found) result.team_layer = std::move(global_solver.best_color);
    } else if (strict_found) {
        result.team_layer = std::move(insertion_solver.best_color);
    }

    if (strict_found) {
        int moved_teams = 0;
        long long moved_nets = 0;
        for (int i = 0; i < n; ++i) {
            if (result.team_layer[i] != baseline_layer[i]) {
                moved_teams++;
                moved_nets += team_size[i];
            }
        }
        std::cout << "[LayerCap] strict " << max_layers << "-layer coloring found, "
                  << (preserved_existing ? "overflow-only insertion" : "global recoloring")
                  << ", search nodes " << strict_nodes
                  << ", moved " << moved_teams << " teams / " << moved_nets << " nets"
                  << std::endl;
        return result;
    }

    // Generic fallback: preserve legal baseline assignments, place overflow teams with
    // minimum conflict, then improve with single moves and pair swaps. This is used only
    // when strict coloring is not found within the bounded search budget.
    result.team_layer.assign(n, -1);
    std::vector<long long> layer_load(max_layers, 0);
    for (int i = 0; i < n; ++i) {
        if (baseline_layer[i] >= 0 && baseline_layer[i] < max_layers) {
            result.team_layer[i] = baseline_layer[i];
            layer_load[baseline_layer[i]] += team_size[i];
        }
    }
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (conflict_graph[a].size() != conflict_graph[b].size()) {
            return conflict_graph[a].size() > conflict_graph[b].size();
        }
        if (team_size[a] != team_size[b]) return team_size[a] > team_size[b];
        return a < b;
    });
    for (int idx : order) {
        if (result.team_layer[idx] >= 0) continue;
        std::tuple<long long, long long, long long, int> best_key{
            LLONG_MAX, LLONG_MAX, LLONG_MAX, INT_MAX};
        int best_layer = 0;
        for (int layer = 0; layer < max_layers; ++layer) {
            long long conflict_edges = 0;
            long long conflict_weight = 0;
            for (int neighbor : conflict_graph[idx]) {
                if (result.team_layer[neighbor] == layer) {
                    conflict_edges++;
                    conflict_weight += team_size[idx] + team_size[neighbor];
                }
            }
            auto key = std::make_tuple(conflict_edges, conflict_weight, layer_load[layer], layer);
            if (key < best_key) {
                best_key = key;
                best_layer = layer;
            }
        }
        result.team_layer[idx] = best_layer;
        layer_load[best_layer] += team_size[idx];
    }

    auto score = [&]() {
        long long conflict_edges = 0;
        long long conflict_weight = 0;
        long long moved_nets = 0;
        std::vector<long long> loads(max_layers, 0);
        for (int i = 0; i < n; ++i) {
            loads[result.team_layer[i]] += team_size[i];
            if (result.team_layer[i] != baseline_layer[i]) moved_nets += team_size[i];
            for (int neighbor : conflict_graph[i]) {
                if (i < neighbor && result.team_layer[i] == result.team_layer[neighbor]) {
                    conflict_edges++;
                    conflict_weight += team_size[i] + team_size[neighbor];
                }
            }
        }
        long long load_square = 0;
        for (long long load : loads) load_square += load * load;
        return std::make_tuple(conflict_edges, conflict_weight, moved_nets, load_square);
    };

    for (int iter = 0; iter < 50; ++iter) {
        bool changed = false;
        for (int idx : order) {
            int old_layer = result.team_layer[idx];
            auto best_score = score();
            int best_layer = old_layer;
            for (int layer = 0; layer < max_layers; ++layer) {
                if (layer == old_layer) continue;
                result.team_layer[idx] = layer;
                auto candidate_score = score();
                if (candidate_score < best_score) {
                    best_score = candidate_score;
                    best_layer = layer;
                }
            }
            result.team_layer[idx] = best_layer;
            changed |= best_layer != old_layer;
        }
        for (int a = 0; a < n; ++a) {
            for (int b = a + 1; b < n; ++b) {
                if (result.team_layer[a] == result.team_layer[b]) continue;
                auto before = score();
                std::swap(result.team_layer[a], result.team_layer[b]);
                if (score() < before) {
                    changed = true;
                } else {
                    std::swap(result.team_layer[a], result.team_layer[b]);
                }
            }
        }
        if (!changed) break;
    }
    auto final_score = score();
    std::cout << "[LayerCap] strict coloring not found within budget (nodes " << strict_nodes
              << "), fallback conflicts " << std::get<0>(final_score)
              << ", conflict weight " << std::get<1>(final_score)
              << ", moved nets " << std::get<2>(final_score) << std::endl;
    return result;
}

LayerAssignmentResult LayerAssigner::assign_layers(std::vector<PairTeam>& vpt, RoutingCase& rc) {
    // Apply the same five-grid expansion later used by ExitBuilder.
    auto extend_bounds = [&](PairTeam& pt) {
        pt.start_left_bound_x = std::max(0, pt.start_left_bound_x - TEAM_BOUNDARY_EXTEND);
        pt.start_right_bound_x = std::min(rc.grid_max_width - 1, pt.start_right_bound_x + TEAM_BOUNDARY_EXTEND);
        pt.start_up_bound_y = std::max(0, pt.start_up_bound_y - TEAM_BOUNDARY_EXTEND);
        pt.start_down_bound_y = std::min(rc.grid_max_height - 1, pt.start_down_bound_y + TEAM_BOUNDARY_EXTEND);
        pt.end_left_bound_x = std::max(0, pt.end_left_bound_x - TEAM_BOUNDARY_EXTEND);
        pt.end_right_bound_x = std::min(rc.grid_max_width - 1, pt.end_right_bound_x + TEAM_BOUNDARY_EXTEND);
        pt.end_up_bound_y = std::max(0, pt.end_up_bound_y - TEAM_BOUNDARY_EXTEND);
        pt.end_down_bound_y = std::min(rc.grid_max_height - 1, pt.end_down_bound_y + TEAM_BOUNDARY_EXTEND);
    };
    int next_source_team_id = 0;
    for (auto it = vpt.begin(); it != vpt.end(); ) {
        if (it->cnt == 0) {
            it = vpt.erase(it);
        } else {
            it->source_team_id = next_source_team_id++;
            it->split_family_id = -1;
            extend_bounds(*it);
            ++it;
        }
    }

    LayerAssignmentResult original_result;
    if (vpt.size() <= 20) {
        original_result = assign_optimal_layers(vpt);
    } else {
        original_result = assign_approximate_layers(vpt);
    }
    int preserved_layer_floor = original_result.total_layers;
    if (preserved_layer_floor < 2) preserved_layer_floor = 2;
    if (preserved_layer_floor == 2) preserved_layer_floor = 4;
    preserved_layer_floor = std::min(
        preserved_layer_floor, MAX_ROUTING_LAYERS);

    const uint64_t split_min_work = get_team_split_min_work();
    const uint64_t intersect_split_min_work =
        get_intersect_team_split_min_work();
    const int max_accepted_team_splits = get_max_accepted_team_splits();
    const int max_accepted_intersect_splits =
        get_max_accepted_intersect_team_splits();
    const bool allow_disjoint_sibling_same_layer =
        allow_disjoint_split_siblings_same_layer();
    const bool preserve_parent_layer_for_one_child =
        preserve_parent_split_layer();
    constexpr size_t MIN_TEAM_SPLIT_NETS = 16;
    const size_t original_team_count = vpt.size();
    size_t eligible_teams = 0;
    size_t accepted_splits = 0;
    size_t accepted_intersect_splits = 0;
    size_t accepted_ordinary_splits = 0;
    size_t geometry_rejections = 0;
    size_t layer_rejections = 0;
    LayerAssignmentResult result;
    bool split_result_valid = false;

    std::vector<int> split_order;
    split_order.reserve(vpt.size());
    for (const PairTeam& team : vpt) {
        const bool intersects = endpoint_bounds_overlap(team);
        const uint64_t required_work = intersects
            ? intersect_split_min_work : split_min_work;
        // TEAM_SPLIT_MIN_WORK=0 remains the global off switch used by the
        // paired no-Team-split baseline, irrespective of the intersect knobs.
        if (split_min_work > 0 && max_accepted_team_splits > 0 &&
            required_work > 0 && team.member.size() >= MIN_TEAM_SPLIT_NETS &&
            team_work_proxy(team) >= required_work) {
            split_order.push_back(team.source_team_id);
        }
    }
    std::sort(split_order.begin(), split_order.end(), [&](int lhs, int rhs) {
        const PairTeam& a = vpt[lhs];
        const PairTeam& b = vpt[rhs];
        const bool a_intersects = endpoint_bounds_overlap(a);
        const bool b_intersects = endpoint_bounds_overlap(b);
        if (a_intersects != b_intersects) return a_intersects > b_intersects;
        const uint64_t a_work = team_work_proxy(a);
        const uint64_t b_work = team_work_proxy(b);
        if (a_work != b_work) {
            return a_work > b_work;
        }
        return lhs < rhs;
    });

    std::cout << "[TeamSplitConfig] min_work " << split_min_work
              << ", intersect_min_work " << intersect_split_min_work
              << ", min_nets " << MIN_TEAM_SPLIT_NETS
              << ", max_ordinary_accepted " << max_accepted_team_splits
              << ", max_intersect_accepted "
              << max_accepted_intersect_splits
              << ", sibling_layer_policy "
              << (allow_disjoint_sibling_same_layer
                      ? "disjoint_same"
                      : (preserve_parent_layer_for_one_child
                             ? "preserve_parent" : "different"))
              << ", max_layers " << MAX_ROUTING_LAYERS << std::endl;
    for (int source_team_id : split_order) {
        auto current_it = std::find_if(
            vpt.begin(), vpt.end(), [&](const PairTeam& team) {
                return team.source_team_id == source_team_id &&
                       team.split_family_id < 0;
            });
        if (current_it == vpt.end()) continue;
        const size_t current_idx = static_cast<size_t>(current_it - vpt.begin());
        const PairTeam& parent = *current_it;
        const uint64_t parent_work = team_work_proxy(parent);
        const bool parent_intersects = endpoint_bounds_overlap(parent);
        if (parent_intersects) {
            if (accepted_intersect_splits >=
                static_cast<size_t>(max_accepted_intersect_splits)) {
                continue;
            }
        } else if (accepted_ordinary_splits >=
                   static_cast<size_t>(max_accepted_team_splits)) {
            continue;
        }
        eligible_teams++;
        PairTeam first_child = parent;
        PairTeam second_child = parent;
        if (!build_team_split_candidate(
                parent, rc.grid_max_width, rc.grid_max_height,
                first_child, second_child)) {
            geometry_rejections++;
            std::cout << "[TeamSplit] source " << source_team_id
                      << ", nets " << parent.member.size()
                      << ", work " << parent_work
                      << ", parent_intersects " << parent_intersects
                      << ", rejected geometry" << std::endl;
            continue;
        }
        first_child.source_team_id = source_team_id;
        second_child.source_team_id = source_team_id;
        first_child.split_family_id = source_team_id;
        second_child.split_family_id = source_team_id;
        std::vector<PairTeam> trial = vpt;
        trial[current_idx] = std::move(first_child);
        trial.insert(trial.begin() + current_idx + 1, std::move(second_child));

        const auto trial_graph = build_conflict_graph(trial);
        LayerAssignmentResult trial_result;
        bool preserved_parent_layer = false;
        const LayerAssignmentResult& current_baseline =
            split_result_valid ? result : original_result;
        if (allow_disjoint_sibling_same_layer &&
            current_baseline.team_layer.size() == vpt.size() &&
            current_baseline.total_layers <= MAX_ROUTING_LAYERS) {
            trial_result.team_layer.resize(trial.size(), -1);
            for (size_t trial_idx = 0; trial_idx < trial.size(); ++trial_idx) {
                if (trial_idx < current_idx) {
                    trial_result.team_layer[trial_idx] =
                        current_baseline.team_layer[trial_idx];
                } else if (trial_idx == current_idx ||
                           trial_idx == current_idx + 1) {
                    trial_result.team_layer[trial_idx] =
                        current_baseline.team_layer[current_idx];
                } else {
                    trial_result.team_layer[trial_idx] =
                        current_baseline.team_layer[trial_idx - 1];
                }
            }
            if (count_same_layer_conflicts(
                    trial_graph, trial_result.team_layer) == 0) {
                trial_result.total_layers = std::max(
                    used_layer_count(trial_result.team_layer),
                    preserved_layer_floor);
                trial_result.is_optimal = false;
                preserved_parent_layer = true;
            }
        }
        if (!preserved_parent_layer &&
            preserve_parent_layer_for_one_child &&
            current_baseline.team_layer.size() == vpt.size() &&
            current_baseline.total_layers <= MAX_ROUTING_LAYERS) {
            const int parent_layer = current_baseline.team_layer[current_idx];
            const uint64_t child_work[2] = {
                team_work_proxy(trial[current_idx]),
                team_work_proxy(trial[current_idx + 1])};
            using PreserveScore = std::tuple<uint64_t, int, int, int>;
            PreserveScore best_score{
                std::numeric_limits<uint64_t>::max(), INT_MAX, INT_MAX,
                INT_MAX};
            std::vector<int> best_assignment;
            for (int kept_child = 0; kept_child < 2; ++kept_child) {
                const int moved_child = 1 - kept_child;
                for (int moved_layer = 0;
                     moved_layer < MAX_ROUTING_LAYERS; ++moved_layer) {
                    if (moved_layer == parent_layer) continue;
                    std::vector<int> assignment(trial.size(), -1);
                    for (size_t trial_idx = 0;
                         trial_idx < trial.size(); ++trial_idx) {
                        if (trial_idx < current_idx) {
                            assignment[trial_idx] =
                                current_baseline.team_layer[trial_idx];
                        } else if (trial_idx > current_idx + 1) {
                            assignment[trial_idx] =
                                current_baseline.team_layer[trial_idx - 1];
                        }
                    }
                    assignment[current_idx + kept_child] = parent_layer;
                    assignment[current_idx + moved_child] = moved_layer;
                    if (count_same_layer_conflicts(
                            trial_graph, assignment) != 0) {
                        continue;
                    }
                    const PreserveScore score{
                        child_work[moved_child],
                        std::abs(moved_layer - parent_layer),
                        moved_layer, kept_child};
                    if (best_assignment.empty() || score < best_score) {
                        best_score = score;
                        best_assignment = std::move(assignment);
                    }
                }
            }
            if (!best_assignment.empty()) {
                trial_result.team_layer = std::move(best_assignment);
                trial_result.total_layers = std::max(
                    used_layer_count(trial_result.team_layer),
                    preserved_layer_floor);
                trial_result.is_optimal = false;
                preserved_parent_layer = true;
            }
        }
        if (!preserved_parent_layer) {
            trial_result = assign_approximate_layers(trial);
        }
        if (trial_result.total_layers > MAX_ROUTING_LAYERS) {
            trial_result = assign_bounded_layers(
                trial, trial_graph, trial_result, MAX_ROUTING_LAYERS);
        }
        const int residual_conflicts = count_same_layer_conflicts(
            trial_graph, trial_result.team_layer);
        if (residual_conflicts != 0) {
            layer_rejections++;
            std::cout << "[TeamSplit] source " << source_team_id
                      << ", nets " << parent.member.size()
                      << ", work " << parent_work
                      << ", children " << trial[current_idx].member.size()
                      << "+" << trial[current_idx + 1].member.size()
                      << ", parent_intersects " << parent_intersects
                      << ", rejected 6-layer coloring, conflicts "
                      << residual_conflicts << std::endl;
            continue;
        }

        trial_result.total_layers = std::max(
            used_layer_count(trial_result.team_layer), preserved_layer_floor);
        std::cout << "[TeamSplit] source " << source_team_id
                  << ", nets " << parent.member.size()
                  << ", work " << parent_work
                  << ", children " << trial[current_idx].member.size()
                  << "+" << trial[current_idx + 1].member.size()
                  << ", parent_intersects " << parent_intersects
                  << ", accepted, layers " << trial_result.total_layers
                  << ", sibling_layers "
                  << trial_result.team_layer[current_idx] << "+"
                  << trial_result.team_layer[current_idx + 1]
                  << ", preserved_parent_layer " << preserved_parent_layer
                  << std::endl;
        vpt = std::move(trial);
        result = std::move(trial_result);
        split_result_valid = true;
        accepted_splits++;
        if (parent_intersects) {
            accepted_intersect_splits++;
        } else {
            accepted_ordinary_splits++;
        }
    }

    std::cout << "[TeamSplitSummary] original_teams " << original_team_count
              << ", eligible " << eligible_teams
              << ", accepted " << accepted_splits
              << ", geometry_rejections " << geometry_rejections
              << ", layer_rejections " << layer_rejections
              << ", final_teams " << vpt.size()
              << ", accepted_intersect " << accepted_intersect_splits
              << ", accepted_ordinary " << accepted_ordinary_splits
              << std::endl;

    if (!split_result_valid) {
        result = std::move(original_result);
    }
    if (result.total_layers > MAX_ROUTING_LAYERS) {
        std::vector<int> merge_depth(vpt.size(), 0);
        constexpr int MAX_TEAM_MERGES = 8;

        auto maximum_clique = [](const std::vector<std::vector<int>>& graph) {
            int n = static_cast<int>(graph.size());
            std::vector<std::vector<uint8_t>> adjacent(n, std::vector<uint8_t>(n, 0));
            for (int i = 0; i < n; ++i) {
                for (int neighbor : graph[i]) adjacent[i][neighbor] = 1;
            }
            std::vector<int> best;
            std::vector<int> current;
            std::function<void(std::vector<int>)> expand = [&](std::vector<int> candidates) {
                if (current.size() + candidates.size() <= best.size()) return;
                while (!candidates.empty()) {
                    if (current.size() + candidates.size() <= best.size()) return;
                    int best_pos = 0;
                    int best_degree = -1;
                    for (int pos = 0; pos < static_cast<int>(candidates.size()); ++pos) {
                        int degree = 0;
                        for (int other : candidates) degree += adjacent[candidates[pos]][other] != 0;
                        if (degree > best_degree) {
                            best_degree = degree;
                            best_pos = pos;
                        }
                    }
                    int vertex = candidates[best_pos];
                    candidates.erase(candidates.begin() + best_pos);
                    current.push_back(vertex);
                    std::vector<int> next;
                    next.reserve(candidates.size());
                    for (int candidate : candidates) {
                        if (adjacent[vertex][candidate]) next.push_back(candidate);
                    }
                    if (next.empty()) {
                        if (current.size() > best.size()) best = current;
                    } else {
                        expand(std::move(next));
                    }
                    current.pop_back();
                }
            };
            std::vector<int> candidates(n);
            std::iota(candidates.begin(), candidates.end(), 0);
            expand(std::move(candidates));
            return best;
        };

        auto merge_cost = [&](int a_idx, int b_idx) {
            const PairTeam& a = vpt[a_idx];
            const PairTeam& b = vpt[b_idx];
            int a_left = std::min(a.start_left_bound_x, a.end_left_bound_x);
            int a_right = std::max(a.start_right_bound_x, a.end_right_bound_x);
            int a_up = std::min(a.start_up_bound_y, a.end_up_bound_y);
            int a_down = std::max(a.start_down_bound_y, a.end_down_bound_y);
            int b_left = std::min(b.start_left_bound_x, b.end_left_bound_x);
            int b_right = std::max(b.start_right_bound_x, b.end_right_bound_x);
            int b_up = std::min(b.start_up_bound_y, b.end_up_bound_y);
            int b_down = std::max(b.start_down_bound_y, b.end_down_bound_y);
            long long overlap_width = std::max(0, std::min(a_right, b_right) - std::max(a_left, b_left));
            long long overlap_height = std::max(0, std::min(a_down, b_down) - std::max(a_up, b_up));
            long long overlap_area = overlap_width * overlap_height;
            double cross_direction = std::abs(
                a.aver_unit_vector.first * b.aver_unit_vector.second -
                a.aver_unit_vector.second * b.aver_unit_vector.first);
            cross_direction = std::min(1.0, cross_direction);
            double member_weight = static_cast<double>(a.member.size() + b.member.size());
            double contraction_penalty = 1.0 + 4.0 * (merge_depth[a_idx] + merge_depth[b_idx]);
            return static_cast<double>(overlap_area) * member_weight *
                   (0.2 + cross_direction) * contraction_penalty;
        };

        for (int merge_round = 0; merge_round <= MAX_TEAM_MERGES; ++merge_round) {
            auto conflict_graph = build_conflict_graph(vpt);
            LayerAssignmentResult baseline = assign_approximate_layers(vpt);
            if (baseline.total_layers <= MAX_ROUTING_LAYERS) {
                result = baseline;
                break;
            }
            result = assign_bounded_layers(vpt, conflict_graph, baseline, MAX_ROUTING_LAYERS);
            int residual_conflicts = count_same_layer_conflicts(
                conflict_graph, result.team_layer);
            if (residual_conflicts == 0 || merge_round == MAX_TEAM_MERGES) break;

            std::vector<int> core = maximum_clique(conflict_graph);
            if (core.size() <= static_cast<size_t>(MAX_ROUTING_LAYERS)) {
                core.clear();
                std::vector<uint8_t> in_core(vpt.size(), 0);
                for (int i = 0; i < static_cast<int>(conflict_graph.size()); ++i) {
                    for (int neighbor : conflict_graph[i]) {
                        if (i < neighbor && result.team_layer[i] == result.team_layer[neighbor]) {
                            if (!in_core[i]) {
                                in_core[i] = 1;
                                core.push_back(i);
                            }
                            if (!in_core[neighbor]) {
                                in_core[neighbor] = 1;
                                core.push_back(neighbor);
                            }
                        }
                    }
                }
            }
            int merge_a = -1;
            int merge_b = -1;
            double best_cost = std::numeric_limits<double>::infinity();
            bool found_fresh_pair = false;
            for (size_t i = 0; i < core.size(); ++i) {
                for (size_t j = i + 1; j < core.size(); ++j) {
                    int a = core[i];
                    int b = core[j];
                    bool fresh_pair = merge_depth[a] == 0 && merge_depth[b] == 0;
                    if (found_fresh_pair && !fresh_pair) continue;
                    double cost = merge_cost(a, b);
                    if ((fresh_pair && !found_fresh_pair) || cost < best_cost ||
                        (std::abs(cost - best_cost) <= 1e-6 && std::make_pair(a, b) < std::make_pair(merge_a, merge_b))) {
                        found_fresh_pair = fresh_pair;
                        best_cost = cost;
                        merge_a = a;
                        merge_b = b;
                    }
                }
            }
            if (merge_a < 0 || merge_b < 0) break;
            if (merge_a > merge_b) std::swap(merge_a, merge_b);

            PairTeam merged = vpt[merge_a];
            std::vector<BumpPair> incoming = vpt[merge_b].member;
            double direction_dot =
                vpt[merge_a].aver_unit_vector.first * vpt[merge_b].aver_unit_vector.first +
                vpt[merge_a].aver_unit_vector.second * vpt[merge_b].aver_unit_vector.second;
            if (direction_dot < 0.0) {
                for (BumpPair& bp : incoming) std::swap(bp.first, bp.second);
            }
            size_t first_size = merged.member.size();
            size_t second_size = incoming.size();
            merged.member.insert(merged.member.end(), incoming.begin(), incoming.end());
            rebuild_team_geometry(
                merged, rc.grid_max_width, rc.grid_max_height);
            int new_depth = merge_depth[merge_a] + merge_depth[merge_b] + 1;
            vpt[merge_a] = std::move(merged);
            vpt.erase(vpt.begin() + merge_b);
            merge_depth[merge_a] = new_depth;
            merge_depth.erase(merge_depth.begin() + merge_b);
            std::cout << "[LayerMerge] contracted teams " << merge_a << " and " << merge_b
                      << " (" << first_size << " + " << second_size << " nets), cost "
                      << best_cost << ", previous clique " << core.size() << std::endl;
        }
        auto final_graph = build_conflict_graph(vpt);
        int final_conflicts = count_same_layer_conflicts(
            final_graph, result.team_layer);
        result.residual_same_layer_conflicts = final_conflicts;
        std::cout << "[LayerMerge] final teams " << vpt.size()
                  << ", residual same-layer box conflicts " << final_conflicts << std::endl;
    }
    if (result.total_layers < 2) {
        result.total_layers = 2;
    }
    size_t axis_aligned_two_pin_nets = 0;
    for (const std::string& net_name : rc.net_name_list) {
        const Net& net = rc.net_dict.at(net_name);
        if (net.bump_names.size() != 2) continue;
        const Bump& first = rc.bump_dict.at(net.bump_names[0]);
        const Bump& second = rc.bump_dict.at(net.bump_names[1]);
        const int dx = std::abs(second.grid_coord_x - first.grid_coord_x);
        const int dy = std::abs(second.grid_coord_y - first.grid_coord_y);
        if ((dx == 0) != (dy == 0)) ++axis_aligned_two_pin_nets;
    }
    // Use the full legal budget only when enough axis-aligned nets can
    // amortize the extra layers. Small/non-axis cases keep their original
    // layer count and quality.
    if (axis_aligned_two_pin_nets >= 1000 && result.total_layers < 6) {
        result.total_layers = 6;
    } else if (result.total_layers == 2) {
        result.total_layers = 4;
    }
    // 始终应用层重映射（原来只在 total_layers > 2 时执行）
    {
        rc.layer_count = result.total_layers;
        for (std::string net_name : rc.net_name_list) {
            Net& net = rc.net_dict.at(net_name);
            for (std::string bump_name : net.bump_names) {
                Bump& bump = rc.bump_dict.at(bump_name);
                if (bump.layer == 1) {
                    bump.layer = rc.layer_count - 1;
                }
            }
        }
        for (PairTeam& pt : vpt) {
            for (auto& bump_pair : pt.member) {
                if (bump_pair.first.layer == 1) {
                    bump_pair.first.layer = rc.layer_count - 1;
                }
                if (bump_pair.second.layer == 1) {
                    bump_pair.second.layer = rc.layer_count - 1;
                }
            }
        }
    }
    // 输出到文件时编号从0开始
    const std::string output_dir = rc.output_dir.empty() ? "output" : rc.output_dir;
    const std::string output_filename = output_dir + "/" + rc.case_name + "_layer_assignment.txt";
    std::ofstream outfile(output_filename);
    if (!outfile.is_open()) {
       // std::cerr << "\nError: Failed to open output file " << output_filename << std::endl;
       return result;
    }
    outfile << "PairTeam Layer Assignment Info\n";
    outfile << "---\n";
    for (unsigned int i = 0; i < result.team_layer.size(); i++) {
        outfile << "PairTeamID: " << i << '\n';                  // 改为 i（从0开始）
        outfile << "Layer: " << result.team_layer[i] << '\n';
        outfile << "SourceTeamID: " << vpt[i].source_team_id << '\n';
        outfile << "SplitFamilyID: " << vpt[i].split_family_id << '\n';
        // 添加网络名称输出
        for (const auto& bp : vpt[i].member) {
            outfile << "NetName: " << bp.first.net_name << '\n';
        }
        outfile << "StartLeft: " << vpt[i].start_left_bound_x << '\n';
        outfile << "StartRight: " << vpt[i].start_right_bound_x << '\n';
        outfile << "StartUp: " << vpt[i].start_up_bound_y << '\n';
        outfile << "StartDown: " << vpt[i].start_down_bound_y << '\n';
        outfile << "EndLeft: " << vpt[i].end_left_bound_x << '\n';
        outfile << "EndRight: " << vpt[i].end_right_bound_x << '\n';
        outfile << "EndUp: " << vpt[i].end_up_bound_y << '\n';
        outfile << "EndDown: " << vpt[i].end_down_bound_y << '\n';
        outfile << "---\n";
    }
    outfile.close();
    return result;
}
