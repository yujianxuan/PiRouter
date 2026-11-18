#include "a_star.h"
#include "a_star_state_grid.h"
#include "encode_coord.h"
#include <omp.h>
#include <chrono>
#include <deque>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdlib>
#include <memory_resource>
#include <memory>
#include <numeric>
#include <sstream>

extern std::ofstream log_file;

// ============================================================================
// AStarNode implementation
// ============================================================================

AStarNode::AStarNode(int x_, int y_, int layer_)
    : x(x_), y(y_), layer(layer_), g(0), h(0), f(0), parent(nullptr) {}

bool AStarNode::operator>(const AStarNode& other) const {
    return f > other.f;
}

int grid_distance(int x1, int y1, int x2, int y2) {
    return std::max(abs(x1 - x2), abs(y1 - y2));
}

bool is_point_in_convex_polygon(double x, double y, const std::vector<Point>& edges) {
    int n = edges.size();
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        const Point& p1 = edges[i];
        const Point& p2 = edges[j];
        double cross = (p2.first - p1.first) * (y - p1.second) - (p2.second - p1.second) * (x - p1.first);
        double dist_p1 = sqrt(pow(x - p1.first, 2) + pow(y - p1.second, 2));
        double dist_p2 = sqrt(pow(x - p2.first, 2) + pow(y - p2.second, 2));
        double dist_p1p2 = sqrt(pow(p2.first - p1.first, 2) + pow(p2.second - p1.second, 2));
        if (fabs(cross) < 1e-9 && (dist_p1 + dist_p2 - dist_p1p2) < 1e-9) return true;
        if (((p1.second > y) != (p2.second > y)) &&
            (x < (p2.first - p1.first) * (y - p1.second) / (p2.second - p1.second) + p1.first)) {
            inside = !inside;
        }
    }
    return inside;
}

bool is_same_direction(int dx1, int dy1, int dl1,
                       int dx2, int dy2, int dl2) {
    return (dx1 == dx2) && (dy1 == dy2) && (dl1 == dl2);
}

std::vector<std::pair<ResultPoint, ResultPoint>> split_to_segments(const std::vector<AStarNode>& path) {
    std::vector<std::pair<ResultPoint, ResultPoint>> segments;
    if (path.size() < 2) return segments;

    ResultPoint current_start = {path[0].x, path[0].y, path[0].layer};
    int prev_dx = path[1].x - path[0].x;
    int prev_dy = path[1].y - path[0].y;
    int prev_dl = path[1].layer - path[0].layer;

    for (size_t i = 2; i < path.size(); ++i) {
        const AStarNode& prev_node = path[i-1];
        const AStarNode& curr_node = path[i];
        int curr_dx = curr_node.x - prev_node.x;
        int curr_dy = curr_node.y - prev_node.y;
        int curr_dl = curr_node.layer - prev_node.layer;

        if (!is_same_direction(prev_dx, prev_dy, prev_dl,
                               curr_dx, curr_dy, curr_dl)) {
            ResultPoint new_node = {prev_node.x, prev_node.y, prev_node.layer};
            segments.emplace_back(current_start, new_node);
            current_start = new_node;
            prev_dx = curr_dx; prev_dy = curr_dy; prev_dl = curr_dl;
        }
    }

    const AStarNode& last_node = path.back();
    ResultPoint end_point = {last_node.x, last_node.y, last_node.layer};
    segments.emplace_back(current_start, end_point);
    
    // 验证所有线段方向合法性
    for (const auto& seg : segments) {
        int dx = seg.second.x - seg.first.x;
        int dy = seg.second.y - seg.first.y;
        if (dx != 0 && dy != 0 && std::abs(dx) != std::abs(dy)) {
        }
    }
    return segments;
}

// ============================================================================
// Thread-configurable escape routing (per-team local copy)
// ============================================================================

using LocalSet = std::pmr::unordered_set<uint64_t>;

struct TeamEscapeResult {
    std::unordered_map<std::string, std::vector<Seg>> escape_results;
    std::vector<std::pair<size_t, int>> fail_bumps;
    size_t pattern_routes = 0;
    size_t limited_searches = 0;
    double elapsed_s = 0.0;
};

enum class EscapeWorkKind {
    WHOLE_TEAM,
    START_SIDE,
    END_SIDE
};

struct EscapeWorkItem {
    size_t team_idx;
    EscapeWorkKind kind;
    uint64_t estimated_cost;
    size_t active_pairs;
    size_t pair_begin;
    size_t pair_end;
};

static uint64_t escape_scheduling_cost(const EscapeWorkItem& item) {
    // Geometric distance alone underestimates wide buses whose few residual
    // non-pattern routes can trigger expensive searches.  Add a small floor
    // for the original bus width, but cap it by the remaining active work so
    // a nearly-completed very wide bus cannot masquerade as the dominant
    // tail. Search limits and deterministic commit semantics are unchanged.
    static constexpr uint64_t PER_PAIR_SCHEDULING_FLOOR = 100;
    static constexpr size_t ACTIVE_PAIR_WIDTH_CAP = 40;
    const size_t total_pairs = item.pair_end - item.pair_begin;
    const size_t scheduling_pairs = std::min(
        total_pairs, item.active_pairs * ACTIVE_PAIR_WIDTH_CAP);
    return item.estimated_cost +
        PER_PAIR_SCHEDULING_FLOOR * scheduling_pairs;
}

struct ThreadEscapeLoad {
    std::vector<size_t> work_item_indices;
    size_t net_count = 0;
    size_t escape_task_count = 0;
    size_t whole_team_count = 0;
    size_t half_team_count = 0;
    double worker_elapsed_s = 0.0;
};

static const char* work_kind_name(EscapeWorkKind kind) {
    switch (kind) {
        case EscapeWorkKind::WHOLE_TEAM: return "whole";
        case EscapeWorkKind::START_SIDE: return "A";
        case EscapeWorkKind::END_SIDE: return "B";
    }
    return "unknown";
}

static std::string format_work_items(
    const std::vector<size_t>& work_item_indices,
    const std::vector<EscapeWorkItem>& work_items) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < work_item_indices.size(); ++i) {
        if (i != 0) oss << ",";
        const EscapeWorkItem& item = work_items[work_item_indices[i]];
        oss << item.team_idx << work_kind_name(item.kind);
    }
    oss << "]";
    return oss.str();
}

struct AxisCandidate {
    std::string net_name;
    const Bump* first;
    const Bump* second;
    int length;
    int assigned_layer;
};

struct AxisRouteChoice {
    int layer = -1;
    std::vector<std::pair<int, int>> corners;
};

static bool axis_route_choice_is_clear(
    const AxisCandidate& candidate,
    const AxisRouteChoice& choice,
    const OccGrid& occupied,
    const std::unordered_set<uint64_t>& direct_path_occupied) {

    if (choice.layer < 0 || choice.corners.size() < 2) return false;
    for (const Bump* endpoint : {candidate.first, candidate.second}) {
        const int low = std::min(endpoint->layer, choice.layer);
        const int high = std::max(endpoint->layer, choice.layer);
        for (int layer = low; layer <= high; ++layer) {
            if (direct_path_occupied.count(encode_coord(
                    endpoint->grid_coord_x,
                    endpoint->grid_coord_y,
                    layer))) {
                return false;
            }
        }
    }

    const Bump& second = *candidate.second;
    for (size_t segment_index = 1;
         segment_index < choice.corners.size(); ++segment_index) {
        const int x = choice.corners[segment_index - 1].first;
        const int y = choice.corners[segment_index - 1].second;
        const int end_x = choice.corners[segment_index].first;
        const int end_y = choice.corners[segment_index].second;
        const bool ends_at_final_endpoint =
            end_x == second.grid_coord_x &&
            end_y == second.grid_coord_y;
        if (!occupied.planar_segment_clear(
                x, y, end_x, end_y, choice.layer,
                ends_at_final_endpoint)) {
            return false;
        }
    }
    return true;
}

static AxisRouteChoice find_axis_route_choice(
    const AxisCandidate& candidate,
    const RoutingCase& rc,
    const OccGrid& occupied,
    const std::unordered_set<uint64_t>& direct_path_occupied) {

    const Bump& first = *candidate.first;
    const Bump& second = *candidate.second;
    std::vector<int> layers(rc.layer_count);
    std::iota(layers.begin(), layers.end(), 0);
    std::stable_sort(layers.begin(), layers.end(), [&](int a, int b) {
        const int a_vias = std::abs(first.layer - a) +
                           std::abs(second.layer - a);
        const int b_vias = std::abs(first.layer - b) +
                           std::abs(second.layer - b);
        if (a_vias != b_vias) return a_vias < b_vias;
        const bool a_assigned = a == candidate.assigned_layer;
        const bool b_assigned = b == candidate.assigned_layer;
        if (a_assigned != b_assigned) return a_assigned;
        return a > b;
    });

    for (int layer : layers) {
        bool endpoint_vias_clear = true;
        for (const Bump* endpoint : {&first, &second}) {
            const int low = std::min(endpoint->layer, layer);
            const int high = std::max(endpoint->layer, layer);
            for (int via_layer = low; via_layer <= high; ++via_layer) {
                if (direct_path_occupied.count(encode_coord(
                        endpoint->grid_coord_x,
                        endpoint->grid_coord_y,
                        via_layer))) {
                    endpoint_vias_clear = false;
                    break;
                }
            }
            if (!endpoint_vias_clear) break;
        }
        if (!endpoint_vias_clear) continue;

        const bool axis_aligned =
            first.grid_coord_x == second.grid_coord_x ||
            first.grid_coord_y == second.grid_coord_y;
        static constexpr int RECTILINEAR_OFFSET_LIMIT = 256;
        static constexpr int DIAGONAL_OFFSET_LIMIT = 64;
        const int rectilinear_pattern_count =
            4 * RECTILINEAR_OFFSET_LIMIT;
        const int pattern_count = axis_aligned
            ? 65
            : 2 + rectilinear_pattern_count +
                4 * DIAGONAL_OFFSET_LIMIT;
        for (int pattern_index = 0;
             pattern_index < pattern_count; ++pattern_index) {
            std::vector<std::pair<int, int>> corners;
            corners.emplace_back(first.grid_coord_x, first.grid_coord_y);
            if (axis_aligned) {
                const int offset = pattern_index == 0 ? 0 :
                    ((pattern_index + 1) / 2) *
                    (pattern_index % 2 ? 1 : -1);
                if (first.grid_coord_y == second.grid_coord_y) {
                    const int shifted_y = first.grid_coord_y + offset;
                    if (shifted_y < 0 ||
                        shifted_y >= rc.grid_max_height) continue;
                    if (offset != 0) {
                        corners.emplace_back(
                            first.grid_coord_x, shifted_y);
                        corners.emplace_back(
                            second.grid_coord_x, shifted_y);
                    }
                } else {
                    const int shifted_x = first.grid_coord_x + offset;
                    if (shifted_x < 0 ||
                        shifted_x >= rc.grid_max_width) continue;
                    if (offset != 0) {
                        corners.emplace_back(
                            shifted_x, first.grid_coord_y);
                        corners.emplace_back(
                            shifted_x, second.grid_coord_y);
                    }
                }
            } else if (pattern_index < 2) {
                const int dx = second.grid_coord_x - first.grid_coord_x;
                const int dy = second.grid_coord_y - first.grid_coord_y;
                const int diagonal = std::min(std::abs(dx), std::abs(dy));
                const int sx = (dx > 0) - (dx < 0);
                const int sy = (dy > 0) - (dy < 0);
                if (pattern_index == 0) {
                    corners.emplace_back(
                        first.grid_coord_x + sx * diagonal,
                        first.grid_coord_y + sy * diagonal);
                } else {
                    corners.emplace_back(
                        second.grid_coord_x - sx * diagonal,
                        second.grid_coord_y - sy * diagonal);
                }
            } else if (pattern_index <
                           2 + rectilinear_pattern_count) {
                const int offset_index = pattern_index - 2;
                const int mode = offset_index % 4;
                const int offset = (offset_index / 4 + 1) *
                    (mode % 2 == 0 ? 1 : -1);
                if (std::abs(second.grid_coord_x - first.grid_coord_x) >=
                    std::abs(second.grid_coord_y - first.grid_coord_y)) {
                    const int anchor_y = mode < 2
                        ? first.grid_coord_y : second.grid_coord_y;
                    const int shifted_y = anchor_y + offset;
                    if (shifted_y < 0 ||
                        shifted_y >= rc.grid_max_height) continue;
                    corners.emplace_back(first.grid_coord_x, shifted_y);
                    corners.emplace_back(second.grid_coord_x, shifted_y);
                } else {
                    const int anchor_x = mode < 2
                        ? first.grid_coord_x : second.grid_coord_x;
                    const int shifted_x = anchor_x + offset;
                    if (shifted_x < 0 ||
                        shifted_x >= rc.grid_max_width) continue;
                    corners.emplace_back(shifted_x, first.grid_coord_y);
                    corners.emplace_back(shifted_x, second.grid_coord_y);
                }
            } else if (std::abs(
                           second.grid_coord_x - first.grid_coord_x) >=
                       std::abs(
                           second.grid_coord_y - first.grid_coord_y)) {
                const int offset_index = pattern_index - 2 -
                    rectilinear_pattern_count;
                const int mode = offset_index % 4;
                const int offset = (offset_index / 4 + 1) *
                    (mode % 2 == 0 ? 1 : -1);
                const int anchor_y = mode < 2
                    ? first.grid_coord_y : second.grid_coord_y;
                const int shifted_y = anchor_y + offset;
                if (shifted_y < 0 ||
                    shifted_y >= rc.grid_max_height) continue;
                const int sx =
                    (second.grid_coord_x > first.grid_coord_x) -
                    (second.grid_coord_x < first.grid_coord_x);
                const int first_diagonal =
                    std::abs(shifted_y - first.grid_coord_y);
                const int second_diagonal =
                    std::abs(second.grid_coord_y - shifted_y);
                const int first_trunk_x =
                    first.grid_coord_x + sx * first_diagonal;
                const int second_trunk_x =
                    second.grid_coord_x - sx * second_diagonal;
                if ((second_trunk_x - first_trunk_x) * sx < 0) continue;
                corners.emplace_back(first_trunk_x, shifted_y);
                corners.emplace_back(second_trunk_x, shifted_y);
            } else {
                const int offset_index = pattern_index - 2 -
                    rectilinear_pattern_count;
                const int mode = offset_index % 4;
                const int offset = (offset_index / 4 + 1) *
                    (mode % 2 == 0 ? 1 : -1);
                const int anchor_x = mode < 2
                    ? first.grid_coord_x : second.grid_coord_x;
                const int shifted_x = anchor_x + offset;
                if (shifted_x < 0 ||
                    shifted_x >= rc.grid_max_width) continue;
                const int sy =
                    (second.grid_coord_y > first.grid_coord_y) -
                    (second.grid_coord_y < first.grid_coord_y);
                const int first_diagonal =
                    std::abs(shifted_x - first.grid_coord_x);
                const int second_diagonal =
                    std::abs(second.grid_coord_x - shifted_x);
                const int first_trunk_y =
                    first.grid_coord_y + sy * first_diagonal;
                const int second_trunk_y =
                    second.grid_coord_y - sy * second_diagonal;
                if ((second_trunk_y - first_trunk_y) * sy < 0) continue;
                corners.emplace_back(shifted_x, first_trunk_y);
                corners.emplace_back(shifted_x, second_trunk_y);
            }
            corners.emplace_back(second.grid_coord_x, second.grid_coord_y);
            corners.erase(
                std::unique(corners.begin(), corners.end()), corners.end());

            AxisRouteChoice choice;
            choice.layer = layer;
            choice.corners = std::move(corners);
            if (axis_route_choice_is_clear(
                    candidate, choice, occupied,
                    direct_path_occupied)) {
                return choice;
            }
        }
    }
    return AxisRouteChoice{};
}

// ----------------------------------------------------------------------------
// Neighbor generation with local set and global copy
// ----------------------------------------------------------------------------

struct NeighborCoord {
    int x;
    int y;
    int layer;
};

static size_t get_neighbors_local(
    const AStarNode& current,
    int max_width, int max_height,
    int min_width, int min_height,
    const OccGrid& global_copy,
    const LocalSet& local_occupied,
    int end_x, int end_y, int end_layer,
    std::array<NeighborCoord, 8>& neighbors) {
    size_t neighbor_count = 0;
    static constexpr int dirs[8][2] = {
        {0,1},{0,-1},{1,0},{-1,0},{1,1},{1,-1},{-1,1},{-1,-1}
    };
    for (const auto& dir : dirs) {
        int dx = dir[0], dy = dir[1];
        int new_x = current.x + dx;
        int new_y = current.y + dy;
        if (new_x >= min_width && new_x <= max_width && new_y >= min_height && new_y <= max_height) {
            uint64_t coord = encode_coord(new_x, new_y, end_layer);
            if (global_copy.find(coord) != global_copy.end() ||
                local_occupied.find(coord) != local_occupied.end()) {
                if (!(new_x == end_x && new_y == end_y)) continue;
            }
            if (dx != 0 && dy != 0) {
                if (current.parent != nullptr) {
                    if ((new_x == current.parent->x && current.y == current.parent->y) ||
                        (new_y == current.parent->y && current.x == current.parent->x))
                        continue;
                }
                float mid_x = current.x + dx * 0.5f;
                float mid_y = current.y + dy * 0.5f;
                uint64_t key = encode_coord(mid_x, mid_y, end_layer);
                if (global_copy.find(key) != global_copy.end() ||
                    local_occupied.find(key) != local_occupied.end())
                    continue;
            } else {
                if (dx != 0 && current.parent != nullptr && new_x == current.parent->x) continue;
                if (dy != 0 && current.parent != nullptr && new_y == current.parent->y) continue;
            }
            neighbors[neighbor_count++] = NeighborCoord{new_x, new_y, end_layer};
        }
    }
    return neighbor_count;
}

// ----------------------------------------------------------------------------
// Intersect version of neighbor generation
// ----------------------------------------------------------------------------

static size_t get_neighbors_intersect_local(
    const AStarNode& current,
    const std::vector<Point>& edges,
    const OccGrid& global_copy,
    const LocalSet& local_occupied,
    int end_x, int end_y, int end_layer,
    std::array<NeighborCoord, 8>& neighbors) {
    size_t neighbor_count = 0;
    static constexpr int dirs[8][2] = {
        {0,1},{0,-1},{1,0},{-1,0},{1,1},{1,-1},{-1,1},{-1,-1}
    };
    for (const auto& dir : dirs) {
        int dx = dir[0], dy = dir[1];
        int new_x = current.x + dx;
        int new_y = current.y + dy;
        if (is_point_in_convex_polygon(new_x, new_y, edges)) {
            uint64_t coord = encode_coord(new_x, new_y, end_layer);
            if (global_copy.find(coord) != global_copy.end() ||
                local_occupied.find(coord) != local_occupied.end()) {
                if (!(new_x == end_x && new_y == end_y)) continue;
            }
            if (dx != 0 && dy != 0) {
                if (current.parent != nullptr) {
                    if ((new_x == current.parent->x && current.y == current.parent->y) ||
                        (new_y == current.parent->y && current.x == current.parent->x))
                        continue;
                }
                float mid_x = current.x + dx * 0.5f;
                float mid_y = current.y + dy * 0.5f;
                uint64_t key = encode_coord(mid_x, mid_y, end_layer);
                if (global_copy.find(key) != global_copy.end() ||
                    local_occupied.find(key) != local_occupied.end())
                    continue;
            } else {
                if (dx != 0 && current.parent != nullptr && new_x == current.parent->x) continue;
                if (dy != 0 && current.parent != nullptr && new_y == current.parent->y) continue;
            }
            neighbors[neighbor_count++] = NeighborCoord{new_x, new_y, end_layer};
        }
    }
    return neighbor_count;
}

// ----------------------------------------------------------------------------
// Path reconstruction (inserts into local set)
// ----------------------------------------------------------------------------

static std::vector<AStarNode> reconstruct_path_local(
    AStarNode* end_node,
    LocalSet& local_occupied) {
    
    std::vector<AStarNode> path;
    AStarNode* current = end_node;
    AStarNode* last = end_node;
    while (current) {
        path.push_back(*current);
        uint64_t key = encode_coord(current->x, current->y, current->layer);
        local_occupied.insert(key);
        if (last->x != current->x && last->y != current->y) {
            float mid_x = (last->x + current->x) * 0.5f;
            float mid_y = (last->y + current->y) * 0.5f;
            local_occupied.insert(encode_coord(mid_x, mid_y, current->layer));
        }
        last = current;
        current = current->parent;
    }
    return path;
}

// Try short octilinear paths and nearby rectilinear tracks before invoking
// unit-grid A*.
// This is especially effective for ordered buses: the escape construction has
// already aligned endpoints, so most paths need no obstacle-driven detour.
// The complete path is validated before any occupancy is committed.
template <typename PointAllowed>
static bool try_monotone_escape_pattern(
    const Bump& bump,
    int end_x, int end_y, int end_layer,
    const OccGrid& global_copy,
    LocalSet& local_occupied,
    PointAllowed point_allowed,
    AStarResult& routing_result) {

    const int start_x = bump.grid_coord_x;
    const int start_y = bump.grid_coord_y;
    const int dx_total = end_x - start_x;
    const int dy_total = end_y - start_y;
    const int step_x = (dx_total > 0) - (dx_total < 0);
    const int step_y = (dy_total > 0) - (dy_total < 0);
    const int diagonal_steps = std::min(std::abs(dx_total), std::abs(dy_total));
    const int x_axis_steps = std::abs(dx_total) - diagonal_steps;
    const int y_axis_steps = std::abs(dy_total) - diagonal_steps;

    static constexpr int TRACK_OFFSET_LIMIT = 32;
    const int candidate_count = 2 + 2 * TRACK_OFFSET_LIMIT;
    for (int candidate_index = 0;
         candidate_index < candidate_count; ++candidate_index) {
        std::vector<AStarNode> candidate;
        candidate.reserve(static_cast<size_t>(
            std::abs(dx_total) + std::abs(dy_total) +
            2 * TRACK_OFFSET_LIMIT + 1));
        candidate.emplace_back(start_x, start_y, end_layer);
        int x = start_x;
        int y = start_y;
        auto append_steps = [&](int count, int move_x, int move_y) {
            for (int step = 0; step < count; ++step) {
                x += move_x;
                y += move_y;
                candidate.emplace_back(x, y, end_layer);
            }
        };
        if (candidate_index == 0) {
            append_steps(diagonal_steps, step_x, step_y);
            append_steps(x_axis_steps, step_x, 0);
            append_steps(y_axis_steps, 0, step_y);
        } else if (candidate_index == 1) {
            append_steps(x_axis_steps, step_x, 0);
            append_steps(y_axis_steps, 0, step_y);
            append_steps(diagonal_steps, step_x, step_y);
        } else {
            const int offset_index = candidate_index - 2;
            const int offset_magnitude = offset_index / 2 + 1;
            const int offset = (offset_index % 2 == 0)
                ? offset_magnitude : -offset_magnitude;
            if (std::abs(dx_total) >= std::abs(dy_total)) {
                const int trunk_y = start_y + offset;
                append_steps(std::abs(trunk_y - y), 0,
                             (trunk_y > y) - (trunk_y < y));
                append_steps(std::abs(end_x - x),
                             (end_x > x) - (end_x < x), 0);
                append_steps(std::abs(end_y - y), 0,
                             (end_y > y) - (end_y < y));
            } else {
                const int trunk_x = start_x + offset;
                append_steps(std::abs(trunk_x - x),
                             (trunk_x > x) - (trunk_x < x), 0);
                append_steps(std::abs(end_y - y), 0,
                             (end_y > y) - (end_y < y));
                append_steps(std::abs(end_x - x),
                             (end_x > x) - (end_x < x), 0);
            }
        }

        bool clear = true;
        for (size_t index = 1; index < candidate.size(); ++index) {
            const AStarNode& previous = candidate[index - 1];
            const AStarNode& current = candidate[index];
            if (!point_allowed(current.x, current.y)) {
                clear = false;
                break;
            }
            const bool is_endpoint =
                current.x == end_x && current.y == end_y;
            const uint64_t point_key =
                encode_coord(current.x, current.y, end_layer);
            if (!is_endpoint &&
                (global_copy.count(point_key) ||
                 local_occupied.count(point_key))) {
                clear = false;
                break;
            }
            if (current.x != previous.x && current.y != previous.y) {
                const uint64_t midpoint_key = encode_coord(
                    (current.x + previous.x) * 0.5f,
                    (current.y + previous.y) * 0.5f,
                    end_layer);
                if (global_copy.count(midpoint_key) ||
                    local_occupied.count(midpoint_key)) {
                    clear = false;
                    break;
                }
            }
        }
        if (!clear) continue;

        for (size_t index = 0; index < candidate.size(); ++index) {
            const AStarNode& current = candidate[index];
            local_occupied.insert(
                encode_coord(current.x, current.y, end_layer));
            if (index > 0) {
                const AStarNode& previous = candidate[index - 1];
                if (current.x != previous.x && current.y != previous.y) {
                    local_occupied.insert(encode_coord(
                        (current.x + previous.x) * 0.5f,
                        (current.y + previous.y) * 0.5f,
                        end_layer));
                }
            }
        }
        std::reverse(candidate.begin(), candidate.end());
        routing_result.path = std::move(candidate);
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Local A* search using global copy (no modification to global set)
// ----------------------------------------------------------------------------

static AStarResult a_star_find_path_to_exit_local(
    const Bump& bump,
    int max_width, int max_height,
    int min_width, int min_height,
    const OccGrid& global_copy,
    LocalSet& local_occupied,
    int end_x, int end_y, int end_layer,
    size_t node_limit) {

    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> open_list;
    static thread_local AStarStateGrid search_state;
    search_state.reset(max_width + 1, max_height + 1, end_layer + 1);
    std::deque<AStarNode> parent_nodes;
    AStarResult routing_result;

    AStarNode start_node(bump.grid_coord_x, bump.grid_coord_y, end_layer);
    start_node.h = grid_distance(start_node.x, start_node.y, end_x, end_y);
    start_node.f = start_node.g + 1.1 * start_node.h;
    open_list.push(start_node);
    search_state.set_best(
        start_node.x, start_node.y, start_node.layer, 0);

    size_t nodes_explored = 0;
    while (!open_list.empty()) {
        AStarNode current = open_list.top();
        open_list.pop();

        if (node_limit > 0 && ++nodes_explored > node_limit) {
            routing_result.hit_node_limit = true;
            return routing_result;
        }

        if (current.x == end_x && current.y == end_y && current.layer == end_layer) {
            routing_result.path = reconstruct_path_local(&current, local_occupied);
            return routing_result;
        }

        if (search_state.is_closed(current.x, current.y, current.layer) &&
            search_state.best_g(current.x, current.y, current.layer) < current.g) continue;
        search_state.close(current.x, current.y, current.layer);

        std::array<NeighborCoord, 8> neighbors;
        size_t neighbor_count = get_neighbors_local(
            current, max_width, max_height, min_width, min_height,
            global_copy, local_occupied, end_x, end_y, end_layer, neighbors);
        for (size_t neighbor_index = 0;
             neighbor_index < neighbor_count; ++neighbor_index) {
            const NeighborCoord& coord = neighbors[neighbor_index];
            AStarNode neighbor(coord.x, coord.y, coord.layer);
            if (search_state.is_closed(neighbor.x, neighbor.y, neighbor.layer)) continue;

            int tentative_g = current.g + 1;
            if (search_state.improve(neighbor.x, neighbor.y, neighbor.layer, tentative_g)) {
                parent_nodes.push_back(current);
                neighbor.parent = &parent_nodes.back();
                neighbor.g = tentative_g;
                neighbor.h = grid_distance(neighbor.x, neighbor.y, end_x, end_y);
                neighbor.f = neighbor.g + 1.1 * neighbor.h;
                open_list.push(neighbor);
            }
        }
    }

    routing_result.path = {};
    return routing_result;
}

// ----------------------------------------------------------------------------
// Intersect version of local A* search
// ----------------------------------------------------------------------------

static AStarResult a_star_find_path_to_exit_intersect_local(
    const Bump& bump,
    const std::vector<Point>& edges,
    const OccGrid& global_copy,
    LocalSet& local_occupied,
    int end_x, int end_y, int end_layer,
    size_t node_limit) {

    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> open_list;
    int state_width = std::max(bump.grid_coord_x, end_x) + 1;
    int state_height = std::max(bump.grid_coord_y, end_y) + 1;
    for (const Point& point : edges) {
        state_width = std::max(state_width,
                               static_cast<int>(std::ceil(point.first)) + 1);
        state_height = std::max(state_height,
                                static_cast<int>(std::ceil(point.second)) + 1);
    }
    static thread_local AStarStateGrid search_state;
    search_state.reset(state_width, state_height, end_layer + 1);
    std::deque<AStarNode> parent_nodes;
    AStarResult routing_result;
    // Search radius and node limit

    AStarNode start_node(bump.grid_coord_x, bump.grid_coord_y, end_layer);
    start_node.h = grid_distance(start_node.x, start_node.y, end_x, end_y);
    start_node.f = start_node.g + 1.1 * start_node.h;
    open_list.push(start_node);
    search_state.set_best(
        start_node.x, start_node.y, start_node.layer, 0);

    size_t nodes_explored = 0;
    while (!open_list.empty()) {
        AStarNode current = open_list.top();
        open_list.pop();


        if (node_limit > 0 && ++nodes_explored > node_limit) {
            routing_result.hit_node_limit = true;
            return routing_result;
        }

        if (current.x == end_x && current.y == end_y && current.layer == end_layer) {
            routing_result.path = reconstruct_path_local(&current, local_occupied);
            return routing_result;
        }

        if (search_state.is_closed(current.x, current.y, current.layer) &&
            search_state.best_g(current.x, current.y, current.layer) < current.g) continue;
        search_state.close(current.x, current.y, current.layer);

        std::array<NeighborCoord, 8> neighbors;
        size_t neighbor_count = get_neighbors_intersect_local(
            current, edges, global_copy, local_occupied,
            end_x, end_y, end_layer, neighbors);
        for (size_t neighbor_index = 0;
             neighbor_index < neighbor_count; ++neighbor_index) {
            const NeighborCoord& coord = neighbors[neighbor_index];
            AStarNode neighbor(coord.x, coord.y, coord.layer);
            if (search_state.is_closed(neighbor.x, neighbor.y, neighbor.layer)) continue;

            int tentative_g = current.g + 1;
            if (search_state.improve(neighbor.x, neighbor.y, neighbor.layer, tentative_g)) {
                parent_nodes.push_back(current);
                neighbor.parent = &parent_nodes.back();
                neighbor.g = tentative_g;
                neighbor.h = grid_distance(neighbor.x, neighbor.y, end_x, end_y);
                neighbor.f = neighbor.g + 1.1 * neighbor.h;
                open_list.push(neighbor);
            }
        }
    }

    routing_result.path = {};
    return routing_result;
}

// ----------------------------------------------------------------------------
// Non-overlapping teams can route their two endpoint boxes independently.
// ExitBuilder already computed intersect_flag from these exact bounds, so the
// work queue reuses that decision instead of repeating geometric checks.
// ----------------------------------------------------------------------------

static bool pair_first_is_net_first(
    const BumpPair& bump_pair,
    const RoutingCase& rc) {
    const std::string& net_name = bump_pair.first.net_name;
    const Bump& net_first = rc.bump_dict.at(
        rc.net_dict.at(net_name).bump_names[0]);
    return net_first.grid_coord_x == bump_pair.first.grid_coord_x &&
           net_first.grid_coord_y == bump_pair.first.grid_coord_y &&
           net_first.layer == bump_pair.first.layer;
}

static void a_star_escape_non_intersect_side_local(
    const PairTeam& team,
    bool route_start_side,
    const RoutingCase& rc,
    const std::unordered_set<std::string>& direct_nets,
    LocalSet& local_occupied,
    const OccGrid& global_copy,
    TeamEscapeResult& side_result,
    bool enable_pattern_escape,
    size_t pair_begin,
    size_t pair_end,
    const std::unordered_set<size_t>* skipped_pairs = nullptr) {

    const auto side_start = std::chrono::high_resolution_clock::now();
    const int min_w = route_start_side
        ? team.start_left_bound_x : team.end_left_bound_x;
    const int max_w = route_start_side
        ? team.start_right_bound_x : team.end_right_bound_x;
    const int min_h = route_start_side
        ? team.start_up_bound_y : team.end_up_bound_y;
    const int max_h = route_start_side
        ? team.start_down_bound_y : team.end_down_bound_y;
    static constexpr size_t LARGE_BUS_ASTAR_NODE_LIMIT = 200000;
    const size_t node_limit = enable_pattern_escape
        ? LARGE_BUS_ASTAR_NODE_LIMIT : 0;

    auto endpoint_for_pair = [&](size_t pair_idx) {
        const BumpPair& bump_pair = team.member[pair_idx];
        const std::string& net_name = bump_pair.first.net_name;
        const std::vector<std::string>& bump_names =
            rc.net_dict.at(net_name).bump_names;
        const bool first_is_net_first =
            pair_first_is_net_first(bump_pair, rc);
        const int bump_position = route_start_side
            ? (first_is_net_first ? 0 : 1)
            : (first_is_net_first ? 1 : 0);
        return std::make_pair(
            &rc.bump_dict.at(bump_names[bump_position]), bump_position);
    };

    for (size_t pair_idx = pair_begin; pair_idx < pair_end; ++pair_idx) {
        if (skipped_pairs != nullptr && skipped_pairs->count(pair_idx)) {
            continue;
        }
        const BumpPair& bump_pair = team.member[pair_idx];
        const std::string& net_name = bump_pair.first.net_name;
        if (direct_nets.count(net_name)) continue;
        const std::vector<std::string>& bump_names =
            rc.net_dict.at(net_name).bump_names;
        const int bump_position = endpoint_for_pair(pair_idx).second;
        const std::string& bump_name = bump_names[bump_position];
        const Bump& bump = rc.bump_dict.at(bump_name);

        AStarResult result;
        const bool used_pattern = enable_pattern_escape &&
            try_monotone_escape_pattern(
                bump,
                bump.exit_terminal.grid_coord_x,
                bump.exit_terminal.grid_coord_y,
                bump.exit_terminal.layer,
                global_copy, local_occupied,
                [&](int x, int y) {
                    return x >= min_w && x <= max_w &&
                           y >= min_h && y <= max_h;
                },
                result);
        if (used_pattern) {
            side_result.pattern_routes++;
        } else {
            result = a_star_find_path_to_exit_local(
                bump, max_w, max_h, min_w, min_h,
                global_copy, local_occupied,
                bump.exit_terminal.grid_coord_x,
                bump.exit_terminal.grid_coord_y,
                bump.exit_terminal.layer,
                node_limit);
        }
        if (result.hit_node_limit) side_result.limited_searches++;
        if (result.path.empty() &&
            !(bump.grid_coord_x == bump.exit_terminal.grid_coord_x &&
              bump.grid_coord_y == bump.exit_terminal.grid_coord_y)) {
            side_result.fail_bumps.emplace_back(pair_idx, bump_position);
        }
        std::vector<Seg> new_segments = split_to_segments(result.path);
        auto& destination = side_result.escape_results[bump_name];
        destination.insert(
            destination.end(), new_segments.begin(), new_segments.end());
    }

    const auto side_end = std::chrono::high_resolution_clock::now();
    side_result.elapsed_s =
        std::chrono::duration<double>(side_end - side_start).count();
}

static uint64_t estimate_escape_work_range(
    const PairTeam& team,
    EscapeWorkKind kind,
    const RoutingCase& rc,
    const std::unordered_set<std::string>& direct_nets,
    size_t pair_begin,
    size_t pair_end) {
    uint64_t estimate = 0;
    const bool route_start_side = kind != EscapeWorkKind::END_SIDE;
    for (size_t pair_idx = pair_begin; pair_idx < pair_end; ++pair_idx) {
        const BumpPair& bump_pair = team.member[pair_idx];
        const std::string& net_name = bump_pair.first.net_name;
        if (direct_nets.count(net_name)) continue;
        const std::vector<std::string>& bump_names =
            rc.net_dict.at(net_name).bump_names;
        const bool first_is_net_first = pair_first_is_net_first(bump_pair, rc);
        const int bump_position = route_start_side
            ? (first_is_net_first ? 0 : 1)
            : (first_is_net_first ? 1 : 0);
        const Bump& bump = rc.bump_dict.at(bump_names[bump_position]);
        estimate += static_cast<uint64_t>(grid_distance(
            bump.grid_coord_x, bump.grid_coord_y,
            bump.exit_terminal.grid_coord_x,
            bump.exit_terminal.grid_coord_y)) + 1;
    }
    return estimate;
}

static uint64_t estimate_escape_work(
    const PairTeam& team,
    EscapeWorkKind kind,
    const RoutingCase& rc,
    const std::unordered_set<std::string>& direct_nets) {
    return estimate_escape_work_range(
        team, kind, rc, direct_nets, 0, team.member.size());
}

static size_t count_active_escape_pairs_range(
    const PairTeam& team,
    const std::unordered_set<std::string>& direct_nets,
    size_t pair_begin,
    size_t pair_end) {

    size_t active_pairs = 0;
    for (size_t pair_idx = pair_begin; pair_idx < pair_end; ++pair_idx) {
        if (!direct_nets.count(team.member[pair_idx].first.net_name)) {
            active_pairs++;
        }
    }
    return active_pairs;
}

static uint64_t get_min_split_estimated_cost() {
    // Joint 8/16-thread tuning selected 20k: it preserves routing output while
    // avoiding both tiny-task overhead and the long tails seen above 35k.
    static constexpr uint64_t DEFAULT_MIN_SPLIT_ESTIMATED_COST = 20000;
    const char* raw_value = std::getenv("HALF_TEAM_MIN_SPLIT_COST");
    if (raw_value == nullptr || raw_value[0] == '\0') {
        return DEFAULT_MIN_SPLIT_ESTIMATED_COST;
    }

    errno = 0;
    char* parse_end = nullptr;
    const unsigned long long parsed = std::strtoull(raw_value, &parse_end, 10);
    if (errno != 0 || raw_value[0] == '-' ||
        parse_end == raw_value || *parse_end != '\0') {
        std::cerr << "Invalid HALF_TEAM_MIN_SPLIT_COST: " << raw_value
                  << std::endl;
        std::exit(2);
    }
    return static_cast<uint64_t>(parsed);
}

// ----------------------------------------------------------------------------
// Per-team worker (uses a read-only global snapshot and a team-local set)
// ----------------------------------------------------------------------------

static void a_star_escape_per_team_range_local(
    const PairTeam& team,
    bool intersect_flag,
    RoutingCase& rc,
    const std::unordered_set<std::string>& direct_nets,
    LocalSet& local_occupied,
    const OccGrid& global_copy,
    TeamEscapeResult& team_result,
    size_t pair_begin,
    size_t pair_end,
    bool enable_pattern_escape,
    const std::unordered_set<size_t>* skipped_pairs = nullptr,
    const std::vector<std::pair<size_t, int>>* preserved_failures = nullptr) {

    auto team_start = std::chrono::high_resolution_clock::now();
    static constexpr size_t LARGE_BUS_ASTAR_NODE_LIMIT = 200000;
    const size_t node_limit = enable_pattern_escape
        ? LARGE_BUS_ASTAR_NODE_LIMIT : 0;
    std::vector<std::pair<size_t, int>> fail_bumps;
    if (preserved_failures != nullptr) {
        fail_bumps = *preserved_failures;
    }
    auto& escape_results = team_result.escape_results;

    if (intersect_flag) {
        int a_min_w = team.start_left_bound_x;
        int a_max_w = team.start_right_bound_x;
        int a_min_h = team.start_up_bound_y;
        int a_max_h = team.start_down_bound_y;
        int b_min_w = team.end_left_bound_x;
        int b_max_w = team.end_right_bound_x;
        int b_min_h = team.end_up_bound_y;
        int b_max_h = team.end_down_bound_y;

        std::vector<Point> edges = generate_cut_rect_vertices(
            calculate_occupancy(
                a_min_w, a_max_w, a_min_h, a_max_h,
                b_min_w, b_max_w, b_min_h, b_max_h
            )
        );

        for (size_t pair_idx = pair_begin; pair_idx < pair_end; ++pair_idx) {
            if (skipped_pairs != nullptr && skipped_pairs->count(pair_idx)) {
                continue;
            }
            const BumpPair& bump_pair = team.member[pair_idx];
            const std::string& net_name = bump_pair.first.net_name;
            if (direct_nets.count(net_name)) continue;
            if (rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_x != bump_pair.first.grid_coord_x ||
                rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_y != bump_pair.first.grid_coord_y ||
                rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).layer != bump_pair.first.layer) {
                const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[1];
                const Bump& a_bump = rc.bump_dict.at(a_bump_name);
                AStarResult a_result;
                const bool used_pattern = enable_pattern_escape &&
                    try_monotone_escape_pattern(
                        a_bump,
                        a_bump.exit_terminal.grid_coord_x,
                        a_bump.exit_terminal.grid_coord_y,
                        a_bump.exit_terminal.layer,
                        global_copy, local_occupied,
                        [&](int x, int y) {
                            return is_point_in_convex_polygon(x, y, edges);
                        },
                        a_result);
                if (used_pattern) {
                    team_result.pattern_routes++;
                } else {
                    a_result = a_star_find_path_to_exit_intersect_local(
                        a_bump, edges, global_copy, local_occupied,
                        a_bump.exit_terminal.grid_coord_x,
                        a_bump.exit_terminal.grid_coord_y,
                        a_bump.exit_terminal.layer,
                        node_limit);
                }
                if (a_result.hit_node_limit) team_result.limited_searches++;
                if (a_result.path.empty() && !(a_bump.grid_coord_x == a_bump.exit_terminal.grid_coord_x && a_bump.grid_coord_y == a_bump.exit_terminal.grid_coord_y)) {
                    fail_bumps.emplace_back(pair_idx, 1);
                }
                std::vector<std::pair<ResultPoint, ResultPoint>> new_segments = split_to_segments(a_result.path);
                escape_results[a_bump_name].insert(escape_results[a_bump_name].end(), new_segments.begin(), new_segments.end());
            } else {
                const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[0];
                const Bump& a_bump = rc.bump_dict.at(a_bump_name);
                AStarResult a_result;
                const bool used_pattern = enable_pattern_escape &&
                    try_monotone_escape_pattern(
                        a_bump,
                        a_bump.exit_terminal.grid_coord_x,
                        a_bump.exit_terminal.grid_coord_y,
                        a_bump.exit_terminal.layer,
                        global_copy, local_occupied,
                        [&](int x, int y) {
                            return is_point_in_convex_polygon(x, y, edges);
                        },
                        a_result);
                if (used_pattern) {
                    team_result.pattern_routes++;
                } else {
                    a_result = a_star_find_path_to_exit_intersect_local(
                        a_bump, edges, global_copy, local_occupied,
                        a_bump.exit_terminal.grid_coord_x,
                        a_bump.exit_terminal.grid_coord_y,
                        a_bump.exit_terminal.layer,
                        node_limit);
                }
                if (a_result.hit_node_limit) team_result.limited_searches++;
                if (a_result.path.empty() && !(a_bump.grid_coord_x == a_bump.exit_terminal.grid_coord_x && a_bump.grid_coord_y == a_bump.exit_terminal.grid_coord_y)) {
                    fail_bumps.emplace_back(pair_idx, 0);
                }
                std::vector<std::pair<ResultPoint, ResultPoint>> new_segments = split_to_segments(a_result.path);
                escape_results[a_bump_name].insert(escape_results[a_bump_name].end(), new_segments.begin(), new_segments.end());
            }
        }
    } else {
        int a_min_w = team.start_left_bound_x;
        int a_max_w = team.start_right_bound_x;
        int a_min_h = team.start_up_bound_y;
        int a_max_h = team.start_down_bound_y;
        int b_min_w = team.end_left_bound_x;
        int b_max_w = team.end_right_bound_x;
        int b_min_h = team.end_up_bound_y;
        int b_max_h = team.end_down_bound_y;

        auto route_in_box = [&](const Bump& bump,
                                int min_w, int max_w,
                                int min_h, int max_h) {
            AStarResult result;
            const bool used_pattern = enable_pattern_escape &&
                try_monotone_escape_pattern(
                    bump,
                    bump.exit_terminal.grid_coord_x,
                    bump.exit_terminal.grid_coord_y,
                    bump.exit_terminal.layer,
                    global_copy, local_occupied,
                    [&](int x, int y) {
                        return x >= min_w && x <= max_w &&
                               y >= min_h && y <= max_h;
                    },
                    result);
            if (used_pattern) {
                team_result.pattern_routes++;
            } else {
                result = a_star_find_path_to_exit_local(
                    bump, max_w, max_h, min_w, min_h,
                    global_copy, local_occupied,
                    bump.exit_terminal.grid_coord_x,
                    bump.exit_terminal.grid_coord_y,
                    bump.exit_terminal.layer,
                    node_limit);
            }
            if (result.hit_node_limit) team_result.limited_searches++;
            return result;
        };

        for (size_t pair_idx = pair_begin; pair_idx < pair_end; ++pair_idx) {
            if (skipped_pairs != nullptr && skipped_pairs->count(pair_idx)) {
                continue;
            }
            const BumpPair& bump_pair = team.member[pair_idx];
            const std::string& net_name = bump_pair.first.net_name;
            if (direct_nets.count(net_name)) continue;
            if (rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_x != bump_pair.first.grid_coord_x ||
                rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_y != bump_pair.first.grid_coord_y ||
                rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).layer != bump_pair.first.layer) {
                const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[1];
                const std::string& b_bump_name = rc.net_dict.at(net_name).bump_names[0];
                const Bump& a_bump = rc.bump_dict.at(a_bump_name);
                const Bump& b_bump = rc.bump_dict.at(b_bump_name);

                AStarResult a_result = route_in_box(
                    a_bump, a_min_w, a_max_w, a_min_h, a_max_h);
                if (a_result.path.empty() && !(a_bump.grid_coord_x == a_bump.exit_terminal.grid_coord_x && a_bump.grid_coord_y == a_bump.exit_terminal.grid_coord_y)) {
                    fail_bumps.emplace_back(pair_idx, 1);
                }
                std::vector<std::pair<ResultPoint, ResultPoint>> new_segments_a = split_to_segments(a_result.path);
                escape_results[a_bump_name].insert(escape_results[a_bump_name].end(), new_segments_a.begin(), new_segments_a.end());

                AStarResult b_result = route_in_box(
                    b_bump, b_min_w, b_max_w, b_min_h, b_max_h);
                if (b_result.path.empty() && !(b_bump.grid_coord_x == b_bump.exit_terminal.grid_coord_x && b_bump.grid_coord_y == b_bump.exit_terminal.grid_coord_y)) {
                    fail_bumps.emplace_back(pair_idx, 0);
                }
                std::vector<std::pair<ResultPoint, ResultPoint>> new_segments_b = split_to_segments(b_result.path);
                escape_results[b_bump_name].insert(escape_results[b_bump_name].end(), new_segments_b.begin(), new_segments_b.end());
            } else {
                const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[0];
                const std::string& b_bump_name = rc.net_dict.at(net_name).bump_names[1];
                const Bump& a_bump = rc.bump_dict.at(a_bump_name);
                const Bump& b_bump = rc.bump_dict.at(b_bump_name);

                AStarResult a_result = route_in_box(
                    a_bump, a_min_w, a_max_w, a_min_h, a_max_h);
                if (a_result.path.empty() && !(a_bump.grid_coord_x == a_bump.exit_terminal.grid_coord_x && a_bump.grid_coord_y == a_bump.exit_terminal.grid_coord_y)) {
                    fail_bumps.emplace_back(pair_idx, 0);
                }
                std::vector<std::pair<ResultPoint, ResultPoint>> new_segments_a = split_to_segments(a_result.path);
                escape_results[a_bump_name].insert(escape_results[a_bump_name].end(), new_segments_a.begin(), new_segments_a.end());

                AStarResult b_result = route_in_box(
                    b_bump, b_min_w, b_max_w, b_min_h, b_max_h);
                if (b_result.path.empty() && !(b_bump.grid_coord_x == b_bump.exit_terminal.grid_coord_x && b_bump.grid_coord_y == b_bump.exit_terminal.grid_coord_y)) {
                    fail_bumps.emplace_back(pair_idx, 1);
                }
                std::vector<std::pair<ResultPoint, ResultPoint>> new_segments_b = split_to_segments(b_result.path);
                escape_results[b_bump_name].insert(escape_results[b_bump_name].end(), new_segments_b.begin(), new_segments_b.end());
            }
        }
    }

    auto team_end = std::chrono::high_resolution_clock::now();
    double team_elapsed_ms = std::chrono::duration<double, std::milli>(team_end - team_start).count();
    team_result.fail_bumps = std::move(fail_bumps);
    team_result.elapsed_s = team_elapsed_ms / 1000.0;
}

static void a_star_escape_per_team_local(
    const PairTeam& team,
    bool intersect_flag,
    RoutingCase& rc,
    const std::unordered_set<std::string>& direct_nets,
    LocalSet& local_occupied,
    const OccGrid& global_copy,
    TeamEscapeResult& team_result,
    bool enable_pattern_escape) {
    a_star_escape_per_team_range_local(
        team, intersect_flag, rc, direct_nets, local_occupied, global_copy,
        team_result, 0, team.member.size(), enable_pattern_escape);
}

// ----------------------------------------------------------------------------
// The only escape-routing entry point. Thread count is the sole runtime variable.
// ----------------------------------------------------------------------------

EscapeOutput a_star_escape_with_threads(
    RoutingCase& rc,
    const std::vector<PairTeam>& vpt,
    std::vector<bool> intersect_flags,
    int escape_threads) {

    const auto escape_function_start =
        std::chrono::high_resolution_clock::now();

    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> escape_results;
    std::unordered_map<std::string, std::vector<Seg>> direct_net_paths;
    std::unordered_set<std::string> direct_nets;
    std::unordered_map<int, std::vector<std::pair<size_t, int>>> all_fail_bumps;

    uint64_t largest_team_work = 0;
    uint64_t second_team_work = 0;
    size_t largest_team_nets = 0;
    bool largest_team_intersects = false;
    for (size_t team_idx = 0; team_idx < vpt.size(); ++team_idx) {
        const EscapeWorkKind kind = intersect_flags[team_idx]
            ? EscapeWorkKind::WHOLE_TEAM : EscapeWorkKind::START_SIDE;
        uint64_t work = estimate_escape_work(
            vpt[team_idx], kind, rc, direct_nets);
        if (intersect_flags[team_idx]) work *= 4;
        if (work > largest_team_work) {
            second_team_work = largest_team_work;
            largest_team_work = work;
            largest_team_nets = vpt[team_idx].member.size();
            largest_team_intersects = intersect_flags[team_idx];
        } else if (work > second_team_work) {
            second_team_work = work;
        }
    }
    const bool enable_near_axis_fast_path =
        largest_team_intersects && largest_team_nets >= 1024 &&
        largest_team_work >= 10000000 && second_team_work > 0 &&
        largest_team_work >= 8 * second_team_work;
    const auto escape_profile_end =
        std::chrono::high_resolution_clock::now();

    // 初始化全局 occupied_set（垂直路径 + 终点）
    for (auto& net_name : rc.net_name_list) {
        const Net& net = rc.net_dict.at(net_name);
        for (const auto& bump_name : net.bump_names) {
            const Bump& bump = rc.bump_dict.at(bump_name);
            if (bump.layer <= bump.exit_terminal.layer) {
                for (int l = bump.layer; l <= bump.exit_terminal.layer; ++l) {
                    rc.occupied_set.insert(encode_coord(bump.grid_coord_x, bump.grid_coord_y, l));
                }
                escape_results[bump_name] = { std::make_pair(
                    ResultPoint{bump.grid_coord_x, bump.grid_coord_y, bump.layer},
                    ResultPoint{bump.grid_coord_x, bump.grid_coord_y, bump.exit_terminal.layer}) };
            } else {
                for (int l = bump.exit_terminal.layer; l <= bump.layer; ++l) {
                    rc.occupied_set.insert(encode_coord(bump.grid_coord_x, bump.grid_coord_y, l));
                }
                escape_results[bump_name] = { std::make_pair(
                    ResultPoint{bump.grid_coord_x, bump.grid_coord_y, bump.layer},
                    ResultPoint{bump.grid_coord_x, bump.grid_coord_y, bump.exit_terminal.layer}) };
            }
            rc.occupied_set.insert(encode_coord(bump.exit_terminal.grid_coord_x,
                                                bump.exit_terminal.grid_coord_y,
                                                bump.exit_terminal.layer));
        }
    }
    const auto escape_initial_occupancy_end =
        std::chrono::high_resolution_clock::now();

    std::vector<AxisCandidate> axis_candidates;
    axis_candidates.reserve(rc.net_name_list.size());
    for (const std::string& net_name : rc.net_name_list) {
        const Net& net = rc.net_dict.at(net_name);
        if (net.bump_names.size() != 2) continue;
        const Bump& first = rc.bump_dict.at(net.bump_names[0]);
        const Bump& second = rc.bump_dict.at(net.bump_names[1]);
        const int dx = std::abs(second.grid_coord_x - first.grid_coord_x);
        const int dy = std::abs(second.grid_coord_y - first.grid_coord_y);
        if (dx == 0 && dy == 0) continue;
        if (dx != 0 && dy != 0 &&
            (!enable_near_axis_fast_path ||
             std::min(dx, dy) * 8 > std::max(dx, dy))) {
            continue;
        }
        axis_candidates.push_back(AxisCandidate{
            net_name, &first, &second, dx + dy,
            first.exit_terminal.layer});
    }
    std::stable_sort(axis_candidates.begin(), axis_candidates.end(),
        [](const AxisCandidate& lhs, const AxisCandidate& rhs) {
            if (lhs.length != rhs.length) return lhs.length > rhs.length;
            return lhs.net_name < rhs.net_name;
        });
    static constexpr size_t MIN_AXIS_PATTERN_NETS = 1000;
    if (axis_candidates.size() < MIN_AXIS_PATTERN_NETS) {
        axis_candidates.clear();
    }
    // Reserve every prospective endpoint column before laying any planar
    // pattern. This prevents an earlier pattern from crossing a later via.
    for (const AxisCandidate& candidate : axis_candidates) {
        for (const Bump* endpoint : {candidate.first, candidate.second}) {
            for (int layer = 0; layer < rc.layer_count; ++layer) {
                rc.occupied_set.insert(encode_coord(
                    endpoint->grid_coord_x,
                    endpoint->grid_coord_y,
                    layer));
            }
        }
    }
    long long accepted_axis_length = 0;
    const size_t max_pattern_nets = enable_near_axis_fast_path
        ? 20000 : 10000;
    std::unordered_set<uint64_t> direct_path_occupied;
    std::vector<AxisRouteChoice> prefetched_axis_choices(
        axis_candidates.size());
    const int axis_parallel_threads = escape_threads <= 0
        ? 8 : std::max(1, escape_threads);
    size_t prefetched_axis_reused = 0;
    size_t prefetched_axis_recomputed = 0;
    // Dense axis-dominated designs invalidate many choices within a large
    // snapshot batch. A fixed smaller canonical window keeps enough work for
    // eight workers while reducing stale recomputation. The decision depends
    // only on the candidate population, never on the requested thread count.
    size_t axis_prefetch_batch_size =
        axis_candidates.size() >= 15000 ? 16 : 32;
    if (const char* raw_batch = std::getenv("AXIS_PREFETCH_BATCH_SIZE")) {
        errno = 0;
        char* parse_end = nullptr;
        const unsigned long parsed = std::strtoul(
            raw_batch, &parse_end, 10);
        if (errno != 0 || parse_end == raw_batch || *parse_end != '\0' ||
            parsed == 0 || parsed > 1024) {
            std::cerr << "Invalid AXIS_PREFETCH_BATCH_SIZE: "
                      << raw_batch << std::endl;
            std::exit(2);
        }
        axis_prefetch_batch_size = static_cast<size_t>(parsed);
    }
    for (size_t batch_begin = 0;
         batch_begin < axis_candidates.size();
         batch_begin += axis_prefetch_batch_size) {
        if (direct_nets.size() >= max_pattern_nets) break;
        const size_t batch_end = std::min(
            axis_candidates.size(),
            batch_begin + axis_prefetch_batch_size);
        #pragma omp parallel for schedule(dynamic, 1) \
            num_threads(axis_parallel_threads)
        for (size_t candidate_idx = batch_begin;
             candidate_idx < batch_end; ++candidate_idx) {
            prefetched_axis_choices[candidate_idx] =
                find_axis_route_choice(
                    axis_candidates[candidate_idx], rc,
                    rc.occupied_set, direct_path_occupied);
        }
        for (size_t candidate_idx = batch_begin;
             candidate_idx < batch_end; ++candidate_idx) {
        if (direct_nets.size() >= max_pattern_nets) break;
        const AxisCandidate& candidate = axis_candidates[candidate_idx];
        const Bump& first = *candidate.first;
        const Bump& second = *candidate.second;
        std::vector<int> layers(rc.layer_count);
        std::iota(layers.begin(), layers.end(), 0);
        std::stable_sort(layers.begin(), layers.end(), [&](int a, int b) {
            const int a_vias = std::abs(first.layer - a) +
                               std::abs(second.layer - a);
            const int b_vias = std::abs(first.layer - b) +
                               std::abs(second.layer - b);
            if (a_vias != b_vias) return a_vias < b_vias;
            const bool a_assigned = a == candidate.assigned_layer;
            const bool b_assigned = b == candidate.assigned_layer;
            if (a_assigned != b_assigned) return a_assigned;
            return a > b;
        });
        int selected_layer = -1;
        std::vector<std::pair<int, int>> selected_corners;
        if (axis_route_choice_is_clear(
                candidate, prefetched_axis_choices[candidate_idx],
                rc.occupied_set, direct_path_occupied)) {
            selected_layer = prefetched_axis_choices[candidate_idx].layer;
            selected_corners = std::move(
                prefetched_axis_choices[candidate_idx].corners);
            prefetched_axis_reused++;
        } else {
            prefetched_axis_recomputed++;
        }
        for (int layer : layers) {
            if (selected_layer >= 0) break;
            bool endpoint_vias_clear = true;
            for (const Bump* endpoint : {&first, &second}) {
                const int low = std::min(endpoint->layer, layer);
                const int high = std::max(endpoint->layer, layer);
                for (int via_layer = low; via_layer <= high; ++via_layer) {
                    if (direct_path_occupied.count(encode_coord(
                            endpoint->grid_coord_x,
                            endpoint->grid_coord_y,
                            via_layer))) {
                        endpoint_vias_clear = false;
                        break;
                    }
                }
                if (!endpoint_vias_clear) break;
            }
            if (!endpoint_vias_clear) continue;
            const bool axis_aligned =
                first.grid_coord_x == second.grid_coord_x ||
                first.grid_coord_y == second.grid_coord_y;
            static constexpr int RECTILINEAR_OFFSET_LIMIT = 256;
            static constexpr int DIAGONAL_OFFSET_LIMIT = 64;
            const int rectilinear_pattern_count =
                4 * RECTILINEAR_OFFSET_LIMIT;
            const int pattern_count = axis_aligned
                ? 65
                : 2 + rectilinear_pattern_count +
                    4 * DIAGONAL_OFFSET_LIMIT;
            for (int pattern_index = 0;
                 pattern_index < pattern_count; ++pattern_index) {
                std::vector<std::pair<int, int>> corners;
                corners.emplace_back(first.grid_coord_x, first.grid_coord_y);
                if (axis_aligned) {
                    const int offset = pattern_index == 0 ? 0 :
                        ((pattern_index + 1) / 2) *
                        (pattern_index % 2 ? 1 : -1);
                    if (first.grid_coord_y == second.grid_coord_y) {
                        const int shifted_y = first.grid_coord_y + offset;
                        if (shifted_y < 0 ||
                            shifted_y >= rc.grid_max_height) continue;
                        if (offset != 0) {
                            corners.emplace_back(
                                first.grid_coord_x, shifted_y);
                            corners.emplace_back(
                                second.grid_coord_x, shifted_y);
                        }
                    } else {
                        const int shifted_x = first.grid_coord_x + offset;
                        if (shifted_x < 0 ||
                            shifted_x >= rc.grid_max_width) continue;
                        if (offset != 0) {
                            corners.emplace_back(
                                shifted_x, first.grid_coord_y);
                            corners.emplace_back(
                                shifted_x, second.grid_coord_y);
                        }
                    }
                } else if (pattern_index < 2) {
                    const int dx = second.grid_coord_x - first.grid_coord_x;
                    const int dy = second.grid_coord_y - first.grid_coord_y;
                    const int diagonal = std::min(std::abs(dx), std::abs(dy));
                    const int sx = (dx > 0) - (dx < 0);
                    const int sy = (dy > 0) - (dy < 0);
                    if (pattern_index == 0) {
                        corners.emplace_back(
                            first.grid_coord_x + sx * diagonal,
                            first.grid_coord_y + sy * diagonal);
                    } else {
                        corners.emplace_back(
                            second.grid_coord_x - sx * diagonal,
                            second.grid_coord_y - sy * diagonal);
                    }
                } else if (pattern_index <
                               2 + rectilinear_pattern_count) {
                    const int offset_index = pattern_index - 2;
                    const int mode = offset_index % 4;
                    const int offset = (offset_index / 4 + 1) *
                        (mode % 2 == 0 ? 1 : -1);
                    if (std::abs(
                            second.grid_coord_x - first.grid_coord_x) >=
                        std::abs(
                            second.grid_coord_y - first.grid_coord_y)) {
                        const int anchor_y = mode < 2
                            ? first.grid_coord_y : second.grid_coord_y;
                        const int shifted_y = anchor_y + offset;
                        if (shifted_y < 0 ||
                            shifted_y >= rc.grid_max_height) continue;
                        corners.emplace_back(
                            first.grid_coord_x, shifted_y);
                        corners.emplace_back(
                            second.grid_coord_x, shifted_y);
                    } else {
                        const int anchor_x = mode < 2
                            ? first.grid_coord_x : second.grid_coord_x;
                        const int shifted_x = anchor_x + offset;
                        if (shifted_x < 0 ||
                            shifted_x >= rc.grid_max_width) continue;
                        corners.emplace_back(
                            shifted_x, first.grid_coord_y);
                        corners.emplace_back(
                            shifted_x, second.grid_coord_y);
                    }
                } else if (std::abs(
                               second.grid_coord_x - first.grid_coord_x) >=
                           std::abs(
                               second.grid_coord_y - first.grid_coord_y)) {
                    const int offset_index = pattern_index - 2 -
                        rectilinear_pattern_count;
                    const int mode = offset_index % 4;
                    const int offset = (offset_index / 4 + 1) *
                        (mode % 2 == 0 ? 1 : -1);
                    const int anchor_y = mode < 2
                        ? first.grid_coord_y : second.grid_coord_y;
                    const int shifted_y = anchor_y + offset;
                    if (shifted_y < 0 || shifted_y >= rc.grid_max_height) continue;
                    const int sx = (second.grid_coord_x > first.grid_coord_x) -
                                   (second.grid_coord_x < first.grid_coord_x);
                    const int first_diagonal =
                        std::abs(shifted_y - first.grid_coord_y);
                    const int second_diagonal =
                        std::abs(second.grid_coord_y - shifted_y);
                    const int first_trunk_x =
                        first.grid_coord_x + sx * first_diagonal;
                    const int second_trunk_x =
                        second.grid_coord_x - sx * second_diagonal;
                    if ((second_trunk_x - first_trunk_x) * sx < 0) continue;
                    corners.emplace_back(first_trunk_x, shifted_y);
                    corners.emplace_back(second_trunk_x, shifted_y);
                } else {
                    const int offset_index = pattern_index - 2 -
                        rectilinear_pattern_count;
                    const int mode = offset_index % 4;
                    const int offset = (offset_index / 4 + 1) *
                        (mode % 2 == 0 ? 1 : -1);
                    const int anchor_x = mode < 2
                        ? first.grid_coord_x : second.grid_coord_x;
                    const int shifted_x = anchor_x + offset;
                    if (shifted_x < 0 || shifted_x >= rc.grid_max_width) continue;
                    const int sy = (second.grid_coord_y > first.grid_coord_y) -
                                   (second.grid_coord_y < first.grid_coord_y);
                    const int first_diagonal =
                        std::abs(shifted_x - first.grid_coord_x);
                    const int second_diagonal =
                        std::abs(second.grid_coord_x - shifted_x);
                    const int first_trunk_y =
                        first.grid_coord_y + sy * first_diagonal;
                    const int second_trunk_y =
                        second.grid_coord_y - sy * second_diagonal;
                    if ((second_trunk_y - first_trunk_y) * sy < 0) continue;
                    corners.emplace_back(shifted_x, first_trunk_y);
                    corners.emplace_back(shifted_x, second_trunk_y);
                }
                corners.emplace_back(second.grid_coord_x, second.grid_coord_y);
                corners.erase(
                    std::unique(corners.begin(), corners.end()), corners.end());

                bool clear = true;
                for (size_t segment_index = 1;
                     segment_index < corners.size() && clear; ++segment_index) {
                    const int x = corners[segment_index - 1].first;
                    const int y = corners[segment_index - 1].second;
                    const int end_x = corners[segment_index].first;
                    const int end_y = corners[segment_index].second;
                    const bool ends_at_final_endpoint =
                        end_x == second.grid_coord_x &&
                        end_y == second.grid_coord_y;
                    clear = rc.occupied_set.planar_segment_clear(
                        x, y, end_x, end_y, layer,
                        ends_at_final_endpoint);
                }
                if (clear) {
                    selected_layer = layer;
                    selected_corners = std::move(corners);
                    break;
                }
            }
            if (selected_layer >= 0) break;
        }
        if (selected_layer < 0) continue;

        std::vector<Seg> path;
        if (first.layer != selected_layer) {
            path.push_back({
                ResultPoint{first.grid_coord_x, first.grid_coord_y, first.layer},
                ResultPoint{first.grid_coord_x, first.grid_coord_y, selected_layer}});
        }
        for (size_t corner_index = 1;
             corner_index < selected_corners.size(); ++corner_index) {
            path.push_back({
                ResultPoint{selected_corners[corner_index - 1].first,
                            selected_corners[corner_index - 1].second,
                            selected_layer},
                ResultPoint{selected_corners[corner_index].first,
                            selected_corners[corner_index].second,
                            selected_layer}});
        }
        if (second.layer != selected_layer) {
            path.push_back({
                ResultPoint{second.grid_coord_x, second.grid_coord_y, selected_layer},
                ResultPoint{second.grid_coord_x, second.grid_coord_y, second.layer}});
        }
        // Planar occupancy is already represented by rc.occupied_set.  This
        // auxiliary set is queried only when validating a later endpoint's
        // vertical layer interval, so retain endpoint via columns rather than
        // duplicating every planar grid point in a hash table.
        for (const Bump* endpoint : {&first, &second}) {
            const int low = std::min(endpoint->layer, selected_layer);
            const int high = std::max(endpoint->layer, selected_layer);
            for (int layer = low; layer <= high; ++layer) {
                direct_path_occupied.insert(encode_coord(
                    endpoint->grid_coord_x,
                    endpoint->grid_coord_y,
                    layer));
            }
        }
        for (const Seg& segment : path) {
            int x = segment.first.x;
            int y = segment.first.y;
            int layer = segment.first.layer;
            const int step_x = (segment.second.x > x) -
                               (segment.second.x < x);
            const int step_y = (segment.second.y > y) -
                               (segment.second.y < y);
            const int step_layer = (segment.second.layer > layer) -
                                   (segment.second.layer < layer);
            while (true) {
                const uint64_t key = encode_coord(x, y, layer);
                rc.occupied_set.insert(key);
                if (x == segment.second.x && y == segment.second.y &&
                    layer == segment.second.layer) {
                    break;
                }
                if (step_x != 0 && step_y != 0) {
                    const uint64_t midpoint_key = encode_coord(
                        x + step_x * 0.5f,
                        y + step_y * 0.5f,
                        layer);
                    rc.occupied_set.insert(midpoint_key);
                }
                x += step_x;
                y += step_y;
                layer += step_layer;
            }
        }
        direct_nets.insert(candidate.net_name);
        direct_net_paths.emplace(candidate.net_name, std::move(path));
        accepted_axis_length += candidate.length;
        }
    }
    log_file << "[AxisPrefetch] threads " << axis_parallel_threads
             << ", batch_size " << axis_prefetch_batch_size
             << ", candidates " << axis_candidates.size()
             << ", reused " << prefetched_axis_reused
             << ", recomputed " << prefetched_axis_recomputed
             << std::endl;
    const auto escape_axis_route_end =
        std::chrono::high_resolution_clock::now();

    // Accepted nets no longer need their old bump-to-Team stubs or exits.
    // Rebuild occupancy from paths that will actually be emitted.
    {
        const int max_coord =
            std::max(rc.grid_max_width, rc.grid_max_height) + 10000;
        rc.occupied_set.init(max_coord, rc.layer_count);
        escape_results.clear();
        for (const std::string& net_name : rc.net_name_list) {
            if (direct_nets.count(net_name)) continue;
            const Net& net = rc.net_dict.at(net_name);
            for (const std::string& bump_name : net.bump_names) {
                const Bump& bump = rc.bump_dict.at(bump_name);
                const int low = std::min(bump.layer, bump.exit_terminal.layer);
                const int high = std::max(bump.layer, bump.exit_terminal.layer);
                for (int layer = low; layer <= high; ++layer) {
                    rc.occupied_set.insert(encode_coord(
                        bump.grid_coord_x, bump.grid_coord_y, layer));
                }
                escape_results[bump_name] = {std::make_pair(
                    ResultPoint{bump.grid_coord_x, bump.grid_coord_y, bump.layer},
                    ResultPoint{bump.grid_coord_x, bump.grid_coord_y,
                                bump.exit_terminal.layer})};
                rc.occupied_set.insert(encode_coord(
                    bump.exit_terminal.grid_coord_x,
                    bump.exit_terminal.grid_coord_y,
                    bump.exit_terminal.layer));
            }
        }
        for (const auto& net_path : direct_net_paths) {
            for (const Seg& segment : net_path.second) {
                int x = segment.first.x;
                int y = segment.first.y;
                int layer = segment.first.layer;
                const int step_x = (segment.second.x > x) -
                                   (segment.second.x < x);
                const int step_y = (segment.second.y > y) -
                                   (segment.second.y < y);
                const int step_layer = (segment.second.layer > layer) -
                                       (segment.second.layer < layer);
                while (true) {
                    rc.occupied_set.insert(encode_coord(x, y, layer));
                    if (x == segment.second.x && y == segment.second.y &&
                        layer == segment.second.layer) {
                        break;
                    }
                    if (step_x != 0 && step_y != 0) {
                        rc.occupied_set.insert(encode_coord(
                            x + step_x * 0.5f,
                            y + step_y * 0.5f,
                            layer));
                    }
                    x += step_x;
                    y += step_y;
                    layer += step_layer;
                }
            }
        }
    }
    const auto escape_rebuild_end =
        std::chrono::high_resolution_clock::now();
    log_file << "[AxisFastPath] eligible " << axis_candidates.size()
             << ", accepted " << direct_nets.size()
             << ", direct length " << accepted_axis_length << std::endl;
    std::cout << "[AxisFastPath] eligible " << axis_candidates.size()
              << ", accepted " << direct_nets.size()
              << ", direct length " << accepted_axis_length << std::endl;

    const size_t team_count = vpt.size();
    const uint64_t min_split_estimated_cost =
        get_min_split_estimated_cost();
    static constexpr uint64_t INTERSECT_WORK_FACTOR = 4;
    std::vector<EscapeWorkItem> work_items;
    work_items.reserve(team_count * 2 + 64);
    std::vector<uint8_t> chunked_team(team_count, 0);
    bool enable_pattern_escape = enable_near_axis_fast_path;
    size_t split_team_count = 0;
    for (size_t team_idx = 0; team_idx < team_count; ++team_idx) {
        const size_t active_pairs = count_active_escape_pairs_range(
            vpt[team_idx], direct_nets, 0,
            vpt[team_idx].member.size());
        if (intersect_flags[team_idx]) {
            const EscapeWorkKind kind = EscapeWorkKind::WHOLE_TEAM;
            const uint64_t estimated_cost =
                estimate_escape_work(vpt[team_idx], kind, rc, direct_nets) *
                INTERSECT_WORK_FACTOR;
            work_items.push_back(EscapeWorkItem{
                team_idx, kind, estimated_cost, active_pairs, 0,
                vpt[team_idx].member.size()});
        } else {
            const EscapeWorkKind start_kind = EscapeWorkKind::START_SIDE;
            const EscapeWorkKind end_kind = EscapeWorkKind::END_SIDE;
            const uint64_t start_cost =
                estimate_escape_work(vpt[team_idx], start_kind, rc, direct_nets);
            const uint64_t end_cost =
                estimate_escape_work(vpt[team_idx], end_kind, rc, direct_nets);
            if (start_cost + end_cost >= min_split_estimated_cost) {
                split_team_count++;
                work_items.push_back(EscapeWorkItem{
                    team_idx, start_kind, start_cost, active_pairs, 0,
                    vpt[team_idx].member.size()});
                work_items.push_back(EscapeWorkItem{
                    team_idx, end_kind, end_cost, active_pairs, 0,
                    vpt[team_idx].member.size()});
            } else {
                work_items.push_back(EscapeWorkItem{
                    team_idx, EscapeWorkKind::WHOLE_TEAM,
                    start_cost + end_cost, active_pairs, 0,
                    vpt[team_idx].member.size()});
            }
        }
    }

    // On large inputs, refine only a sparse dominant tail of wide scheduler
    // items. Independent endpoint boxes use two ranges; an intersecting item
    // may use more short contiguous ranges because it otherwise remains one
    // indivisible critical path. If many items are already dominant, the
    // existing queue has enough parallelism and is left unchanged. Commit
    // below preserves original kind/range order and reconciles overlaps.
    //
    // Decomposition is deliberately independent of the requested worker
    // count.  A 1-thread and an N-thread run therefore execute exactly the
    // same logical work items; only scheduling differs.  Four dominant items
    // is the canonical upper bound for a sparse tail, rather than a value
    // derived from the runtime thread count.
    constexpr size_t MIN_CHUNK_DESIGN_NETS = 3000;
    constexpr size_t MIN_CHUNK_ITEM_PAIRS = 64;
    constexpr uint64_t MIN_CHUNK_ESTIMATED_COST = 20000;
    constexpr uint64_t HIGH_PRESSURE_SCHEDULING_COST = 3000000;
    constexpr size_t MAX_SPARSE_DOMINANT_ITEMS = 4;
    if (rc.net_name_list.size() >= MIN_CHUNK_DESIGN_NETS) {
        uint64_t max_scheduling_cost = 0;
        for (const EscapeWorkItem& item : work_items) {
            max_scheduling_cost = std::max(
                max_scheduling_cost, escape_scheduling_cost(item));
        }
        const uint64_t dominant_cost_cutoff =
            max_scheduling_cost >= HIGH_PRESSURE_SCHEDULING_COST
                ? (max_scheduling_cost + 3) / 4
                : (max_scheduling_cost * 2 + 4) / 5;
        size_t dominant_item_count = 0;
        for (const EscapeWorkItem& item : work_items) {
            if (item.active_pairs >= MIN_CHUNK_ITEM_PAIRS &&
                item.estimated_cost >= MIN_CHUNK_ESTIMATED_COST &&
                escape_scheduling_cost(item) >= dominant_cost_cutoff) {
                dominant_item_count++;
            }
        }
        const bool sparse_dominant_tail =
            dominant_item_count <= MAX_SPARSE_DOMINANT_ITEMS;
        std::vector<EscapeWorkItem> split_items;
        split_items.reserve(work_items.size() * 2);
        for (const EscapeWorkItem& item : work_items) {
            const size_t pair_count = item.pair_end - item.pair_begin;
            if (!sparse_dominant_tail ||
                item.active_pairs < MIN_CHUNK_ITEM_PAIRS ||
                item.estimated_cost < MIN_CHUNK_ESTIMATED_COST ||
                escape_scheduling_cost(item) < dominant_cost_cutoff) {
                split_items.push_back(item);
                continue;
            }
            const size_t chunk_count =
                item.kind == EscapeWorkKind::WHOLE_TEAM &&
                intersect_flags[item.team_idx]
                    ? std::min<size_t>(8, (pair_count + 15) / 16)
                    : 2;
            for (size_t chunk_idx = 0;
                 chunk_idx < chunk_count; ++chunk_idx) {
                const size_t range_begin = item.pair_begin +
                    pair_count * chunk_idx / chunk_count;
                const size_t range_end = item.pair_begin +
                    pair_count * (chunk_idx + 1) / chunk_count;
                uint64_t cost = estimate_escape_work_range(
                    vpt[item.team_idx], item.kind, rc, direct_nets,
                    range_begin, range_end);
                if (item.kind == EscapeWorkKind::WHOLE_TEAM &&
                    intersect_flags[item.team_idx]) {
                    cost *= INTERSECT_WORK_FACTOR;
                }
                split_items.push_back(EscapeWorkItem{
                    item.team_idx, item.kind, cost,
                    count_active_escape_pairs_range(
                        vpt[item.team_idx], direct_nets,
                        range_begin, range_end),
                    range_begin, range_end});
            }
            chunked_team[item.team_idx] = 1;
        }
        work_items.swap(split_items);
    }
    std::stable_sort(
        work_items.begin(), work_items.end(),
        [](const EscapeWorkItem& lhs, const EscapeWorkItem& rhs) {
            const uint64_t lhs_scheduling_cost = escape_scheduling_cost(lhs);
            const uint64_t rhs_scheduling_cost = escape_scheduling_cost(rhs);
            if (lhs_scheduling_cost != rhs_scheduling_cost) {
                return lhs_scheduling_cost > rhs_scheduling_cost;
            }
            if (lhs.team_idx != rhs.team_idx) return lhs.team_idx < rhs.team_idx;
            return lhs.pair_begin < rhs.pair_begin;
        });
    std::vector<std::vector<size_t>> team_work_items(team_count);
    for (size_t work_idx = 0; work_idx < work_items.size(); ++work_idx) {
        team_work_items[work_items[work_idx].team_idx].push_back(work_idx);
    }
    for (std::vector<size_t>& indices : team_work_items) {
        std::sort(
            indices.begin(), indices.end(),
            [&work_items](size_t lhs, size_t rhs) {
                if (work_items[lhs].kind != work_items[rhs].kind) {
                    return static_cast<int>(work_items[lhs].kind) <
                           static_cast<int>(work_items[rhs].kind);
                }
                return work_items[lhs].pair_begin <
                       work_items[rhs].pair_begin;
            });
    }
    const size_t work_item_count = work_items.size();
    const int requested_threads = (escape_threads <= 0)
        ? static_cast<int>(team_count)
        : escape_threads;
    const int active_threads = std::max(
        1, std::min<int>(requested_threads, static_cast<int>(work_item_count)));


    // Each work item owns its allocator, occupancy, and result containers.
    std::vector<std::unique_ptr<std::pmr::unsynchronized_pool_resource>> pools(
        work_item_count);
    for (size_t i = 0; i < work_item_count; ++i) {
        pools[i] = std::make_unique<std::pmr::unsynchronized_pool_resource>();
    }

    std::vector<LocalSet> local_sets;
    local_sets.reserve(work_item_count);
    for (size_t i = 0; i < work_item_count; ++i) {
        local_sets.emplace_back(pools[i].get());
        local_sets.back().reserve(
            work_items[i].kind == EscapeWorkKind::WHOLE_TEAM ? 200000 : 100000);
    }
    const auto escape_queue_prepare_end =
        std::chrono::high_resolution_clock::now();
    std::vector<TeamEscapeResult> work_results(work_item_count);
    std::vector<ThreadEscapeLoad> thread_loads(active_threads);

    omp_set_dynamic(0);
    int observed_threads = 1;
    #pragma omp parallel num_threads(active_threads) shared(observed_threads, thread_loads)
    {
        #pragma omp single
        observed_threads = omp_get_num_threads();

        #pragma omp for schedule(dynamic, 1)
        for (size_t work_idx = 0; work_idx < work_item_count; ++work_idx) {
            const int thread_id = omp_get_thread_num();
            const EscapeWorkItem& item = work_items[work_idx];
            const PairTeam& team = vpt[item.team_idx];
            const bool item_pattern_escape = enable_pattern_escape ||
                (chunked_team[item.team_idx] &&
                 intersect_flags[item.team_idx]);
            if (item.kind == EscapeWorkKind::WHOLE_TEAM) {
                a_star_escape_per_team_range_local(
                    team, intersect_flags[item.team_idx], rc,
                    direct_nets,
                    local_sets[work_idx], rc.occupied_set,
                    work_results[work_idx],
                    item.pair_begin, item.pair_end,
                    item_pattern_escape);
            } else {
                a_star_escape_non_intersect_side_local(
                    team, item.kind == EscapeWorkKind::START_SIDE, rc,
                    direct_nets,
                    local_sets[work_idx], rc.occupied_set,
                    work_results[work_idx],
                    item_pattern_escape,
                    item.pair_begin, item.pair_end);
            }
            ThreadEscapeLoad& load = thread_loads[thread_id];
            load.work_item_indices.push_back(work_idx);
            load.net_count += item.pair_end - item.pair_begin;
            load.escape_task_count += item.pair_end - item.pair_begin;
            if (item.kind == EscapeWorkKind::WHOLE_TEAM) {
                load.whole_team_count++;
            } else {
                load.half_team_count++;
            }
            load.worker_elapsed_s += work_results[work_idx].elapsed_s;
        }
    }
    const auto escape_parallel_end =
        std::chrono::high_resolution_clock::now();
    log_file << "[EscapeWorkQueue] teams " << team_count
             << ", split_teams " << split_team_count
             << ", whole_teams " << (team_count - split_team_count)
             << ", work_items " << work_item_count
             << ", min_split_cost " << min_split_estimated_cost
             << ", intersect_cost_factor " << INTERSECT_WORK_FACTOR
             << std::endl;
    std::cout << "[EscapeWorkQueue] teams " << team_count
              << ", split_teams " << split_team_count
              << ", whole_teams " << (team_count - split_team_count)
              << ", work_items " << work_item_count
              << ", min_split_cost " << min_split_estimated_cost
              << ", intersect_cost_factor " << INTERSECT_WORK_FACTOR
              << std::endl;
    log_file << "[EscapeThreads] requested " << escape_threads
             << ", active " << observed_threads << std::endl;
    std::cout << "[EscapeThreads] requested " << escape_threads
              << ", active " << observed_threads << std::endl;
    for (size_t work_idx = 0; work_idx < work_item_count; ++work_idx) {
        const EscapeWorkItem& item = work_items[work_idx];
        log_file << "[EscapeWorkItem] team " << item.team_idx
                 << ", kind " << work_kind_name(item.kind)
                 << ", nets " << item.pair_end - item.pair_begin
                 << ", range " << item.pair_begin << ":" << item.pair_end
                 << ", estimated_cost " << item.estimated_cost
                 << ", scheduling_cost " << escape_scheduling_cost(item)
                 << ", active_nets " << item.active_pairs
                 << ", pattern_routes " << work_results[work_idx].pattern_routes
                 << ", limited_searches " << work_results[work_idx].limited_searches
                 << ", failed_bumps " << work_results[work_idx].fail_bumps.size()
                 << ", completed in " << work_results[work_idx].elapsed_s
                 << " s" << std::endl;
    }
    for (int thread_id = 0; thread_id < observed_threads; ++thread_id) {
        const ThreadEscapeLoad& load = thread_loads[thread_id];
        std::unordered_set<size_t> unique_teams;
        for (size_t work_idx : load.work_item_indices) {
            unique_teams.insert(work_items[work_idx].team_idx);
        }
        log_file << "[EscapeThreadLoad] thread " << thread_id
                 << " teams " << unique_teams.size()
                 << ", work_items " << load.work_item_indices.size()
                 << ", whole_teams " << load.whole_team_count
                 << ", half_teams " << load.half_team_count
                 << ", nets " << load.net_count
                 << ", escape_tasks " << load.escape_task_count
                 << ", worker_time " << load.worker_elapsed_s << " s"
                 << ", task_ids "
                 << format_work_items(load.work_item_indices, work_items)
                 << std::endl;
    }

    // Merge endpoint tasks back into one deterministic candidate per team.
    // Failure order remains pair-major (A then B), matching the old worker.
    std::vector<TeamEscapeResult> team_results(team_count);
    for (size_t team_idx = 0; team_idx < team_count; ++team_idx) {
        const std::vector<size_t>& indices = team_work_items[team_idx];
        if (chunked_team[team_idx]) continue;
        if (indices.size() == 1) {
            team_results[team_idx] = std::move(work_results[indices[0]]);
            continue;
        }

        TeamEscapeResult& merged = team_results[team_idx];
        TeamEscapeResult& start_result = work_results[indices[0]];
        TeamEscapeResult& end_result = work_results[indices[1]];
        merged.elapsed_s = std::max(
            start_result.elapsed_s, end_result.elapsed_s);
        merged.pattern_routes = start_result.pattern_routes +
                                end_result.pattern_routes;
        merged.limited_searches = start_result.limited_searches +
                                  end_result.limited_searches;
        for (auto& kv : start_result.escape_results) {
            merged.escape_results.emplace(kv.first, std::move(kv.second));
        }
        for (auto& kv : end_result.escape_results) {
            merged.escape_results.emplace(kv.first, std::move(kv.second));
        }

        size_t start_fail_idx = 0;
        size_t end_fail_idx = 0;
        for (size_t pair_idx = 0;
             pair_idx < vpt[team_idx].member.size(); ++pair_idx) {
            if (start_fail_idx < start_result.fail_bumps.size() &&
                start_result.fail_bumps[start_fail_idx].first == pair_idx) {
                merged.fail_bumps.push_back(
                    start_result.fail_bumps[start_fail_idx++]);
            }
            if (end_fail_idx < end_result.fail_bumps.size() &&
                end_result.fail_bumps[end_fail_idx].first == pair_idx) {
                merged.fail_bumps.push_back(
                    end_result.fail_bumps[end_fail_idx++]);
            }
        }
    }

    // Transactional commit: parallel workers search against the same read-only
    // snapshot, so candidates from different teams may overlap. Commit candidates
    // deterministically; reroute only a conflicting team against all paths already
    // committed. This keeps the hot A* loop lock-free without relying on box
    // separation for correctness.
    OccGrid committed_escape;
    int max_coord = std::max(rc.grid_max_width, rc.grid_max_height) + 10000;
    committed_escape.init(max_coord, rc.layer_count);
    int serial_reroute_teams = 0;
    size_t candidate_conflict_cells = 0;

    for (size_t i = 0; i < team_count; ++i) {
        if (chunked_team[i]) {
            TeamEscapeResult& merged = team_results[i];
            double max_parallel_elapsed = 0.0;
            double retry_elapsed = 0.0;
            size_t retried_chunks = 0;
            size_t chunk_conflict_cells = 0;
            for (size_t work_idx : team_work_items[i]) {
                EscapeWorkItem& item = work_items[work_idx];
                max_parallel_elapsed = std::max(
                    max_parallel_elapsed, work_results[work_idx].elapsed_s);
                size_t work_conflict_cells = 0;
                for (uint64_t coord : local_sets[work_idx]) {
                    if (committed_escape.count(coord)) work_conflict_cells++;
                }
                if (work_conflict_cells > 0) {
                    retried_chunks++;
                    chunk_conflict_cells += work_conflict_cells;
                    candidate_conflict_cells += work_conflict_cells;
                    // A failed search produced no occupancy. Adding paths from
                    // earlier chunks cannot make that failed route feasible,
                    // so retain its failure record and replay only candidate
                    // successes. This avoids repeating the most expensive
                    // exhausted searches during transactional reconciliation.
                    const std::vector<std::pair<size_t, int>>
                        preserved_failures =
                            work_results[work_idx].fail_bumps;
                    std::unordered_set<size_t> skipped_pairs;
                    skipped_pairs.reserve(preserved_failures.size());
                    for (const auto& failure : preserved_failures) {
                        skipped_pairs.insert(failure.first);
                    }
                    local_sets[work_idx].clear();
                    work_results[work_idx] = TeamEscapeResult{};
                    const bool item_pattern_escape =
                        enable_pattern_escape ||
                        (chunked_team[i] && intersect_flags[i]);
                    if (item.kind == EscapeWorkKind::WHOLE_TEAM) {
                        a_star_escape_per_team_range_local(
                            vpt[i], intersect_flags[i], rc, direct_nets,
                            local_sets[work_idx], rc.occupied_set,
                            work_results[work_idx],
                            item.pair_begin, item.pair_end,
                            item_pattern_escape,
                            &skipped_pairs, &preserved_failures);
                    } else {
                        work_results[work_idx].fail_bumps =
                            preserved_failures;
                        a_star_escape_non_intersect_side_local(
                            vpt[i],
                            item.kind == EscapeWorkKind::START_SIDE,
                            rc, direct_nets,
                            local_sets[work_idx], rc.occupied_set,
                            work_results[work_idx],
                            item_pattern_escape,
                            item.pair_begin, item.pair_end,
                            &skipped_pairs);
                    }
                    retry_elapsed += work_results[work_idx].elapsed_s;
                }
                for (uint64_t coord : local_sets[work_idx]) {
                    rc.occupied_set.insert(coord);
                    committed_escape.insert(coord);
                }
                for (auto& kv : work_results[work_idx].escape_results) {
                    merged.escape_results.emplace(
                        kv.first, std::move(kv.second));
                }
                merged.fail_bumps.insert(
                    merged.fail_bumps.end(),
                    work_results[work_idx].fail_bumps.begin(),
                    work_results[work_idx].fail_bumps.end());
                merged.pattern_routes +=
                    work_results[work_idx].pattern_routes;
                merged.limited_searches +=
                    work_results[work_idx].limited_searches;
            }
            if (retried_chunks > 0) serial_reroute_teams++;
            merged.elapsed_s = max_parallel_elapsed + retry_elapsed;
            std::stable_sort(
                merged.fail_bumps.begin(), merged.fail_bumps.end(),
                [](const auto& lhs, const auto& rhs) {
                    return lhs.first < rhs.first;
                });
            for (auto& kv : merged.escape_results) {
                auto& dst = escape_results[kv.first];
                dst.insert(dst.end(), kv.second.begin(), kv.second.end());
            }
            all_fail_bumps[static_cast<int>(i)] =
                std::move(merged.fail_bumps);
            log_file << "[EscapeChunkCommit] team " << i
                     << ", chunks " << team_work_items[i].size()
                     << ", retried " << retried_chunks
                     << ", conflict_cells " << chunk_conflict_cells
                     << ", effective_elapsed " << merged.elapsed_s
                     << " s" << std::endl;
            log_file << "[Escape Team " << i << "] completed in "
                     << merged.elapsed_s << " s" << std::endl;
            continue;
        }
        size_t team_conflict_cells = 0;
        for (size_t work_idx : team_work_items[i]) {
            for (uint64_t coord : local_sets[work_idx]) {
                if (committed_escape.count(coord)) team_conflict_cells++;
            }
        }
        if (team_conflict_cells > 0) {
            serial_reroute_teams++;
            candidate_conflict_cells += team_conflict_cells;
            double parallel_elapsed = team_results[i].elapsed_s;
            for (size_t work_idx : team_work_items[i]) {
                local_sets[work_idx].clear();
            }
            team_results[i] = TeamEscapeResult{};
            const size_t fallback_work_idx = team_work_items[i].front();
            a_star_escape_per_team_local(
                vpt[i], intersect_flags[i], rc,
                direct_nets,
                local_sets[fallback_work_idx], rc.occupied_set,
                team_results[i], enable_pattern_escape);
            team_results[i].elapsed_s += parallel_elapsed;
            log_file << "[EscapeCommit] team " << i << " candidate overlapped "
                     << team_conflict_cells << " committed cells; serially rerouted"
                     << std::endl;
        }
        for (size_t work_idx : team_work_items[i]) {
            for (uint64_t coord : local_sets[work_idx]) {
                rc.occupied_set.insert(coord);
                committed_escape.insert(coord);
            }
        }
        for (auto& kv : team_results[i].escape_results) {
            auto& dst = escape_results[kv.first];
            dst.insert(dst.end(), kv.second.begin(), kv.second.end());
        }
        all_fail_bumps[static_cast<int>(i)] = std::move(team_results[i].fail_bumps);
        log_file << "[Escape Team " << i << "] completed in "
                 << team_results[i].elapsed_s << " s" << std::endl;
    }
    log_file << "[EscapeCommit] serially rerouted " << serial_reroute_teams
             << " teams, rejected candidate overlap cells " << candidate_conflict_cells
             << std::endl;
    std::cout << "[EscapeCommit] serially rerouted " << serial_reroute_teams
              << " teams, rejected candidate overlap cells " << candidate_conflict_cells
              << std::endl;

    const auto escape_commit_end =
        std::chrono::high_resolution_clock::now();
    const auto phase_seconds = [](const auto& begin, const auto& end) {
        return std::chrono::duration<double>(end - begin).count();
    };
    log_file << "[EscapePhaseTiming] profile "
             << phase_seconds(escape_function_start, escape_profile_end)
             << " s, initial_occupancy "
             << phase_seconds(escape_profile_end,
                              escape_initial_occupancy_end)
             << " s, axis_route "
             << phase_seconds(escape_initial_occupancy_end,
                              escape_axis_route_end)
             << " s, occupancy_rebuild "
             << phase_seconds(escape_axis_route_end, escape_rebuild_end)
             << " s, queue_prepare "
             << phase_seconds(escape_rebuild_end, escape_queue_prepare_end)
             << " s, parallel_search "
             << phase_seconds(escape_queue_prepare_end, escape_parallel_end)
             << " s, log_and_commit "
             << phase_seconds(escape_parallel_end, escape_commit_end)
             << " s" << std::endl;

    EscapeOutput output;
    output.escape_results = escape_results;
    output.direct_net_paths = std::move(direct_net_paths);
    output.all_fail_bumps = all_fail_bumps;
    return output;
}
