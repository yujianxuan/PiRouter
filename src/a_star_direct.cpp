#include "a_star_direct.h"
#include "a_star_state_grid.h"
#include "encode_coord.h"
#include <deque>
#include "builder.h"
#include "elements.h"
#include <omp.h>
#include <mutex>
#include <shared_mutex>
#include <chrono>
#include <fstream>
#include <cmath>
#include <array>

extern std::ofstream log_file;
extern std::mutex log_mutex;

AStarNodeDirect::AStarNodeDirect(int x_, int y_, int layer_)
    : x(x_), y(y_), layer(layer_), g(0), h(0), f(0), parent(nullptr) {}

bool AStarNodeDirect::operator>(const AStarNodeDirect& other) const {
    return f > other.f;
}

// 删除 get_coord_node_key_direct 和 get_node_key_direct 函数

int grid_distance_direct(int x1, int y1, int l1, int x2, int y2, int l2) {
    return std::max(abs(x1 - x2), abs(y1 - y2)) + abs(l1 - l2);
}

//====================================Single Thread====================================

struct DirectNeighborCoord {
    int x;
    int y;
    int layer;
};

static size_t get_neighbors_direct(
    const AStarNodeDirect& current,
    int max_width,
    int max_height,
    int layer_count,
    const OccGrid& occupied_set,
    int end_x,
    int end_y,
    int end_layer,
    std::array<DirectNeighborCoord, 10>& neighbors) {
    size_t neighbor_count = 0;
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
                if (!(new_x == end_x && new_y == end_y && new_layer == end_layer)) continue;
            }
            if (dx != 0 && dy != 0) {
                if (current.parent != nullptr) {
                    if ((new_x == current.parent->x && current.y == current.parent->y) ||
                        (new_y == current.parent->y && current.x == current.parent->x))
                        continue;
                }
                float mid_x = static_cast<float>(current.x) + static_cast<float>(dx) * 0.5f;
                float mid_y = static_cast<float>(current.y) + static_cast<float>(dy) * 0.5f;
                uint64_t key = encode_coord(mid_x, mid_y, new_layer);
                if (occupied_set.find(key) != occupied_set.end()) continue;
            } else {
                if (dx != 0) {
                    if (current.parent != nullptr && new_x == current.parent->x) continue;
                } else if (dy != 0) {
                    if (current.parent != nullptr && new_y == current.parent->y) continue;
                }
            }
            neighbors[neighbor_count++] =
                DirectNeighborCoord{new_x, new_y, new_layer};
        }
    }
    return neighbor_count;
}

std::vector<AStarNodeDirect> reconstruct_path_direct(
    AStarNodeDirect* end_node,
    OccGrid& occupied_set,
    bool commit_path) {
    std::vector<AStarNodeDirect> path;
    AStarNodeDirect* current = end_node;
    AStarNodeDirect* last = end_node;
    while (current != nullptr) {
        path.push_back(*current);
        uint64_t key = encode_coord(current->x, current->y, current->layer);
        if (commit_path) occupied_set.insert(key);
        if ((last->x != current->x) && (last->y != current->y)) {
            float mid_x = static_cast<float>(last->x + current->x) * 0.5f;
            float mid_y = static_cast<float>(last->y + current->y) * 0.5f;
            uint64_t mid_key = encode_coord(mid_x, mid_y, current->layer);
            if (commit_path) occupied_set.insert(mid_key);
        }
        last = current;
        current = current->parent;
    }
    return path;
}

AStarResultDirect a_star_find_path_direct(const Bump& start,
                                          int max_width, int max_height, int layer_count,
                                          OccGrid& occupied_set,
                                          bool commit_path) {
    const ExitTerminal& target_exit = start.exit_terminal;
    int end_x = target_exit.grid_coord_x;
    int end_y = target_exit.grid_coord_y;
    int end_layer = target_exit.layer;

    std::priority_queue<AStarNodeDirect, std::vector<AStarNodeDirect>, std::greater<AStarNodeDirect>> open_list;
    static thread_local AStarStateGrid search_state;
    search_state.reset(max_width, max_height, layer_count);
    std::deque<AStarNodeDirect> parent_nodes;
    AStarResultDirect routing_result;
    // node limit
    const int MAX_NODES = 2000000;
    int nodes_explored = 0;

    AStarNodeDirect start_node(start.grid_coord_x, start.grid_coord_y, start.layer);
    start_node.h = grid_distance_direct(start.grid_coord_x, start.grid_coord_y, start.layer,
                                        end_x, end_y, end_layer);
    start_node.f = start_node.g + 1.1*start_node.h;
    open_list.push(start_node);
    search_state.set_best(
        start_node.x, start_node.y, start_node.layer, 0);

    while (!open_list.empty()) {
        AStarNodeDirect current = open_list.top();
        open_list.pop();
    nodes_explored++;
    if (nodes_explored > MAX_NODES) {
        routing_result.path = {};
        return routing_result;
    }

        if (current.x == end_x && current.y == end_y && current.layer == end_layer) {
            routing_result.path = reconstruct_path_direct(
                &current, occupied_set, commit_path);
                return routing_result;
        }

        if (search_state.is_closed(current.x, current.y, current.layer) &&
            search_state.best_g(current.x, current.y, current.layer) < current.g) continue;
        search_state.close(current.x, current.y, current.layer);

        std::array<DirectNeighborCoord, 10> neighbors;
        const size_t neighbor_count = get_neighbors_direct(
            current, max_width, max_height, layer_count,
            occupied_set, end_x, end_y, end_layer, neighbors);
        for (size_t neighbor_idx = 0;
             neighbor_idx < neighbor_count; ++neighbor_idx) {
            const DirectNeighborCoord& coord = neighbors[neighbor_idx];
            AStarNodeDirect neighbor(coord.x, coord.y, coord.layer);
            if (search_state.is_closed(neighbor.x, neighbor.y, neighbor.layer)) continue;

            int tentative_g = current.g + 1;
            if (search_state.improve(neighbor.x, neighbor.y, neighbor.layer, tentative_g)) {
                parent_nodes.push_back(current);
                neighbor.parent = &parent_nodes.back();
                neighbor.g = tentative_g;
                neighbor.h = grid_distance_direct(neighbor.x, neighbor.y, neighbor.layer,
                                                  end_x, end_y, end_layer);
                neighbor.f = neighbor.g + 1.1*neighbor.h;
                open_list.push(neighbor);
            }
        }
    }

    routing_result.path = {};
    return routing_result;
}

// 兜底连接专用邻居：与 get_neighbors_direct 相同，但额外允许穿行 self_cells（本 net 自占格点）
// 搜索限制在 [min_x,max_x]×[min_y,max_y] 包围盒内（防止在超大网格上漫游爆 MAX_NODES）
static std::vector<AStarNodeDirect> get_neighbors_connect(const AStarNodeDirect& current,
        int min_x, int max_x, int min_y, int max_y, int layer_count,
        const OccGrid& occupied_set,
        const std::unordered_set<uint64_t>& self_cells,
        int end_x, int end_y, int end_layer) {
    std::vector<AStarNodeDirect> neighbors;
    neighbors.reserve(10);
    static constexpr int dirs[10][3] = {
        {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0},
        {1, 1, 0}, {1, -1, 0}, {-1, 1, 0}, {-1, -1, 0},
        {0, 0, 1}, {0, 0, -1}
    };
    for (const auto& dir : dirs) {
        int dx = dir[0], dy = dir[1], dl = dir[2];
        int new_x = current.x + dx, new_y = current.y + dy, new_layer = current.layer + dl;
        if (new_x >= min_x && new_x <= max_x && new_y >= min_y && new_y <= max_y && new_layer >= 0 && new_layer < layer_count) {
            uint64_t coord = encode_coord(new_x, new_y, new_layer);
            bool is_end = (new_x == end_x && new_y == end_y && new_layer == end_layer);
            if (occupied_set.find(coord) != occupied_set.end()) {
                // 占用格点：仅当是终点、或是本 net 自占格点时才允许穿行
                if (!is_end && self_cells.find(coord) == self_cells.end()) continue;
            }
            if (dx != 0 && dy != 0) {
                if (current.parent != nullptr) {
                    if ((new_x == current.parent->x && current.y == current.parent->y) ||
                        (new_y == current.parent->y && current.x == current.parent->x)) continue;
                }
                float mid_x = static_cast<float>(current.x) + static_cast<float>(dx) * 0.5f;
                float mid_y = static_cast<float>(current.y) + static_cast<float>(dy) * 0.5f;
                uint64_t key = encode_coord(mid_x, mid_y, new_layer);
                if (occupied_set.find(key) != occupied_set.end() && self_cells.find(key) == self_cells.end()) continue;
            } else {
                if (dx != 0) { if (current.parent != nullptr && new_x == current.parent->x) continue; }
                else if (dy != 0) { if (current.parent != nullptr && new_y == current.parent->y) continue; }
            }
            neighbors.emplace_back(new_x, new_y, new_layer);
        }
    }
    return neighbors;
}

AStarResultDirect a_star_connect_points(int sx, int sy, int sl,
                                        int ex, int ey, int el,
                                        int max_width, int max_height, int layer_count,
                                        OccGrid& occupied_set,
                                        const std::unordered_set<uint64_t>& self_cells) {
    int end_x = ex, end_y = ey, end_layer = el;
    // 限制搜索在 start/end 包围盒 + margin 内，避免在超大网格上漫游爆 MAX_NODES
    int margin = std::max(2000, std::max(std::abs(ex - sx), std::abs(ey - sy)));
    int bxmin = std::max(0, std::min(sx, ex) - margin);
    int bxmax = std::min(max_width - 1, std::max(sx, ex) + margin);
    int bymin = std::max(0, std::min(sy, ey) - margin);
    int bymax = std::min(max_height - 1, std::max(sy, ey) + margin);
    std::priority_queue<AStarNodeDirect, std::vector<AStarNodeDirect>, std::greater<AStarNodeDirect>> open_list;
    static thread_local AStarStateGrid search_state;
    search_state.reset(max_width, max_height, layer_count);
    std::deque<AStarNodeDirect> parent_nodes;
    AStarResultDirect routing_result;
    const int MAX_NODES = 30000000;  // 兜底专用：节点预算大（仅极少数失败 bump 调用，不影响 reroute）
    int nodes_explored = 0;

    AStarNodeDirect start_node(sx, sy, sl);
    start_node.h = grid_distance_direct(sx, sy, sl, end_x, end_y, end_layer);
    start_node.f = start_node.g + 1.1*start_node.h;
    open_list.push(start_node);
    search_state.set_best(
        start_node.x, start_node.y, start_node.layer, 0);

    while (!open_list.empty()) {
        AStarNodeDirect current = open_list.top();
        open_list.pop();
        nodes_explored++;
        if (nodes_explored > MAX_NODES) {
            routing_result.path = {};
                return routing_result;
        }

        if (current.x == end_x && current.y == end_y && current.layer == end_layer) {
            routing_result.path = reconstruct_path_direct(&current, occupied_set);
                return routing_result;
        }

        if (search_state.is_closed(current.x, current.y, current.layer) &&
            search_state.best_g(current.x, current.y, current.layer) < current.g) continue;
        search_state.close(current.x, current.y, current.layer);

        std::vector<AStarNodeDirect> neighbors = get_neighbors_connect(current, bxmin, bxmax, bymin, bymax, layer_count,
                                                                      occupied_set, self_cells, end_x, end_y, end_layer);
        for (AStarNodeDirect& neighbor : neighbors) {
            if (search_state.is_closed(neighbor.x, neighbor.y, neighbor.layer)) continue;

            int tentative_g = current.g + 1;
            if (search_state.improve(neighbor.x, neighbor.y, neighbor.layer, tentative_g)) {
                parent_nodes.push_back(current);
                neighbor.parent = &parent_nodes.back();
                neighbor.g = tentative_g;
                neighbor.h = grid_distance_direct(neighbor.x, neighbor.y, neighbor.layer,
                                                  end_x, end_y, end_layer);
                neighbor.f = neighbor.g + 1.1*neighbor.h;
                open_list.push(neighbor);
            }
        }
    }

    routing_result.path = {};
    return routing_result;
}

json coord_to_json_direct(int x, int y, int layer) {
    std::string layer_name = (layer == 0) ? "Bottom" : "Top";
    return json::array({x, y, layer_name});
}

bool is_same_direction_direct(int dx1, int dy1, int dl1,
                              int dx2, int dy2, int dl2) {
    return (dx1 == dx2) && (dy1 == dy2) && (dl1 == dl2);
}

std::vector<std::pair<ResultPoint, ResultPoint>> split_to_segments_direct(const std::vector<AStarNodeDirect>& path) {
    std::vector<std::pair<ResultPoint, ResultPoint>> segments;
    if (path.size() < 2) return segments;

    ResultPoint current_start = {path[0].x, path[0].y, path[0].layer};
    int prev_dx = path[1].x - path[0].x;
    int prev_dy = path[1].y - path[0].y;
    int prev_dl = path[1].layer - path[0].layer;

    for (size_t i = 2; i < path.size(); ++i) {
        const AStarNodeDirect& prev_node = path[i-1];
        const AStarNodeDirect& curr_node = path[i];
        int curr_dx = curr_node.x - prev_node.x;
        int curr_dy = curr_node.y - prev_node.y;
        int curr_dl = curr_node.layer - prev_node.layer;

        if (!is_same_direction_direct(prev_dx, prev_dy, prev_dl,
                                      curr_dx, curr_dy, curr_dl)) {
            ResultPoint new_node = {prev_node.x, prev_node.y, prev_node.layer};
            segments.emplace_back(current_start, new_node);
            current_start = new_node;
            prev_dx = curr_dx; prev_dy = curr_dy; prev_dl = curr_dl;
        }
    }

    const AStarNodeDirect& last_node = path.back();
    ResultPoint end_point = {last_node.x, last_node.y, last_node.layer};
    segments.emplace_back(current_start, end_point);
    return segments;
}

void printFailedNets(std::vector<std::basic_string<char>> &reroute_fail_nets) {
    if (reroute_fail_nets.empty()) return;
    std::cout << "Reroute failed for net: ";
    for (const auto& net_name : reroute_fail_nets) {
        std::cout << net_name << ", ";
    }
    std::cout << "\n";
}

RerouteResult a_star_reroute(
    RoutingCase& rc,
    const std::vector<PairTeam>& vpt,
    std::unordered_map<int, std::vector<std::pair<size_t, int>>> all_fail_bumps) {
    std::vector<std::string> reroute_fail_nets;
    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> escape_results;
    for (size_t team_idx = 0; team_idx < vpt.size(); ++team_idx) {
        auto team_start = std::chrono::high_resolution_clock::now();
        auto fail_iter = all_fail_bumps.find(static_cast<int>(team_idx));
        if (fail_iter == all_fail_bumps.end()) continue;
        const PairTeam& pt = vpt[team_idx];
        for (const auto& fail_pair : fail_iter->second) {
            const BumpPair& bump_pair = pt.member[fail_pair.first];
            const std::string& net_name = bump_pair.first.net_name;
            const std::string& fail_bump_name =
                rc.net_dict.at(net_name).bump_names[fail_pair.second];
            const Bump& fail_bump = rc.bump_dict.at(fail_bump_name);
            AStarResultDirect routing_result = a_star_find_path_direct(
                fail_bump, rc.grid_max_width, rc.grid_max_height,
                rc.layer_count, rc.occupied_set);
            if (routing_result.path.empty() &&
                !(fail_bump.grid_coord_x ==
                      fail_bump.exit_terminal.grid_coord_x &&
                  fail_bump.grid_coord_y ==
                      fail_bump.exit_terminal.grid_coord_y)) {
                reroute_fail_nets.push_back(net_name);
            }
            std::vector<std::pair<ResultPoint, ResultPoint>> new_segments =
                split_to_segments_direct(routing_result.path);
            auto& destination = escape_results[fail_bump_name];
            destination.insert(
                destination.end(), new_segments.begin(), new_segments.end());
        }
        const auto team_end = std::chrono::high_resolution_clock::now();
        const double team_elapsed_s = std::chrono::duration<double>(
            team_end - team_start).count();
        log_file << "[Reroute Team " << team_idx << "] completed in "
                 << team_elapsed_s << " s" << std::endl;
    }
    RerouteResult result;
    result.reroute_fail_nets = reroute_fail_nets;
    result.escape_results = escape_results;
    return result;
}

//====================================Multi Thread====================================

std::vector<AStarNodeDirect> get_neighbors_direct_lock(
    const AStarNodeDirect& current,
    int max_width, int max_height, int layer_count,
    const OccGrid& occupied_set,
    int end_x, int end_y, int end_layer,
    std::shared_mutex& rc_mutex) {
    std::vector<AStarNodeDirect> neighbors;
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
            bool coord_occupied = false;
            {
                std::shared_lock<std::shared_mutex> lock(rc_mutex);
                coord_occupied = (occupied_set.find(coord) != occupied_set.end());
            }
            if (coord_occupied && !(new_x == end_x && new_y == end_y && new_layer == end_layer)) continue;

            if (dx != 0 && dy != 0) {
                if (current.parent != nullptr && ((new_x == current.parent->x && current.y == current.parent->y) ||
                                                  (new_y == current.parent->y && current.x == current.parent->x)))
                    continue;
                float mid_x = static_cast<float>(current.x) + static_cast<float>(dx) * 0.5f;
                float mid_y = static_cast<float>(current.y) + static_cast<float>(dy) * 0.5f;
                uint64_t key = encode_coord(mid_x, mid_y, new_layer);
                bool key_occupied = false;
                {
                    std::shared_lock<std::shared_mutex> lock(rc_mutex);
                    key_occupied = (occupied_set.find(key) != occupied_set.end());
                }
                if (key_occupied) continue;
            } else {
                if (dx != 0) {
                    if (current.parent != nullptr && new_x == current.parent->x) continue;
                } else if (dy != 0) {
                    if (current.parent != nullptr && new_y == current.parent->y) continue;
                }
            }
            neighbors.emplace_back(new_x, new_y, new_layer);
        }
    }
    return neighbors;
}

std::vector<AStarNodeDirect> reconstruct_path_direct_lock(
    AStarNodeDirect* end_node,
    OccGrid& occupied_set,
    std::shared_mutex& rc_mutex) {
    std::vector<AStarNodeDirect> path;
    AStarNodeDirect* current = end_node;
    AStarNodeDirect* last = end_node;
    while (current != nullptr) {
        path.push_back(*current);
        uint64_t key = encode_coord(current->x, current->y, current->layer);
        {
            std::lock_guard<std::shared_mutex> lock(rc_mutex);
            occupied_set.insert(key);
        }
        if ((last->x != current->x) && (last->y != current->y)) {
            float mid_x = static_cast<float>(last->x + current->x) * 0.5f;
            float mid_y = static_cast<float>(last->y + current->y) * 0.5f;
            uint64_t mid_key = encode_coord(mid_x, mid_y, current->layer);
            {
                std::lock_guard<std::shared_mutex> lock(rc_mutex);
                occupied_set.insert(mid_key);
            }
        }
        last = current;
        current = current->parent;
    }
    return path;
}

AStarResultDirect a_star_find_path_direct_lock(
    const Bump& start,
    int max_width, int max_height, int layer_count,
    OccGrid& occupied_set,
    std::shared_mutex& rc_mutex) {
    const ExitTerminal& target_exit = start.exit_terminal;
    int end_x = target_exit.grid_coord_x;
    int end_y = target_exit.grid_coord_y;
    int end_layer = target_exit.layer;

    std::priority_queue<AStarNodeDirect, std::vector<AStarNodeDirect>, std::greater<AStarNodeDirect>> open_list;
    static thread_local AStarStateGrid search_state;
    search_state.reset(max_width, max_height, layer_count);
    std::deque<AStarNodeDirect> parent_nodes;
    AStarResultDirect routing_result;
    // node limit
    const int MAX_NODES = 2000000;
    int nodes_explored = 0;

    AStarNodeDirect start_node(start.grid_coord_x, start.grid_coord_y, start.layer);
    start_node.h = grid_distance_direct(start.grid_coord_x, start.grid_coord_y, start.layer,
                                        end_x, end_y, end_layer);
    start_node.f = start_node.g + 1.1*start_node.h;
    open_list.push(start_node);
    search_state.set_best(
        start_node.x, start_node.y, start_node.layer, 0);

    while (!open_list.empty()) {
        AStarNodeDirect current = open_list.top();
        open_list.pop();
    nodes_explored++;
    if (nodes_explored > MAX_NODES) {
        routing_result.path = {};
        return routing_result;
    }

        if (current.x == end_x && current.y == end_y && current.layer == end_layer) {
            routing_result.path = reconstruct_path_direct_lock(&current, occupied_set, rc_mutex);
                return routing_result;
        }

        if (search_state.is_closed(current.x, current.y, current.layer) &&
            search_state.best_g(current.x, current.y, current.layer) < current.g) continue;
        search_state.close(current.x, current.y, current.layer);

        std::vector<AStarNodeDirect> neighbors = get_neighbors_direct_lock(
            current, max_width, max_height, layer_count, occupied_set,
            end_x, end_y, end_layer, rc_mutex);
        for (AStarNodeDirect& neighbor : neighbors) {
            if (search_state.is_closed(neighbor.x, neighbor.y, neighbor.layer)) continue;

            int tentative_g = current.g + 1;
            if (search_state.improve(neighbor.x, neighbor.y, neighbor.layer, tentative_g)) {
                parent_nodes.push_back(current);
                neighbor.parent = &parent_nodes.back();
                neighbor.g = tentative_g;
                neighbor.h = grid_distance_direct(neighbor.x, neighbor.y, neighbor.layer,
                                                  end_x, end_y, end_layer);
                neighbor.f = neighbor.g + 1.1*neighbor.h;
                open_list.push(neighbor);
            }
        }
    }

    routing_result.path = {};
    return routing_result;
}

void a_star_reroute_per_team(
    RoutingCase& rc,
    const PairTeam& pt,
    const std::vector<std::pair<size_t, int>>& fail_bumps,
    std::vector<std::string> &reroute_fail_nets,
    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> &escape_results,
    std::mutex& fail_nets_mutex,
    std::mutex& escape_results_mutex,
    std::shared_mutex& rc_mutex,
    size_t team_idx) {
    auto team_start = std::chrono::high_resolution_clock::now();

    for (const auto& fail_pair : fail_bumps) {
        size_t pair_idx = fail_pair.first;
        int bump_idx = fail_pair.second;
        const BumpPair& bump_pair = pt.member[pair_idx];
        const std::string& net_name = bump_pair.first.net_name;
        const std::string& fail_bump_name = rc.net_dict.at(net_name).bump_names[bump_idx];
        const Bump& fail_bump = rc.bump_dict.at(fail_bump_name);

        AStarResultDirect routing_result = a_star_find_path_direct_lock(
            fail_bump,
            rc.grid_max_width,
            rc.grid_max_height,
            rc.layer_count,
            rc.occupied_set,
            rc_mutex
        );

        if (routing_result.path.empty() && !(fail_bump.grid_coord_x == fail_bump.exit_terminal.grid_coord_x &&
                                             fail_bump.grid_coord_y == fail_bump.exit_terminal.grid_coord_y)) {
            std::lock_guard<std::mutex> lock(fail_nets_mutex);
            reroute_fail_nets.push_back(net_name);
        }

        std::vector<std::pair<ResultPoint, ResultPoint>> new_segments = split_to_segments_direct(routing_result.path);

        {
            std::lock_guard<std::mutex> lock(escape_results_mutex);
            escape_results[fail_bump_name].insert(escape_results[fail_bump_name].end(),
                                                  new_segments.begin(), new_segments.end());
        }
    }

    auto team_end = std::chrono::high_resolution_clock::now();
    double team_elapsed_ms = std::chrono::duration<double, std::milli>(team_end - team_start).count();
    double team_elapsed_s = team_elapsed_ms / 1000.0;
    {
        std::lock_guard<std::mutex> lock(log_mutex);
        log_file << "[Reroute Team " << team_idx << "] completed in " << team_elapsed_s << " s" << std::endl;
    }
}

RerouteResult a_star_reroute_multithread(
    RoutingCase& rc,
    const std::vector<PairTeam>& vpt,
    std::unordered_map<int, std::vector<std::pair<size_t, int>>> all_fail_bumps) {
    std::vector<std::string> reroute_fail_nets;
    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> escape_results;

    std::mutex fail_nets_mutex;
    std::mutex escape_results_mutex;
    std::shared_mutex rc_mutex;

    omp_set_num_threads(vpt.size());

    size_t team_idx;
    decltype(all_fail_bumps)::iterator fail_bumps_iter;
    const std::vector<std::pair<size_t, int>>* fail_bumps;
    const PairTeam* pt;

    #pragma omp parallel for private(team_idx, fail_bumps_iter, fail_bumps, pt) shared(vpt, all_fail_bumps, rc, reroute_fail_nets, escape_results, fail_nets_mutex, escape_results_mutex, rc_mutex)
    for (team_idx = 0; team_idx < vpt.size(); ++team_idx) {
        fail_bumps_iter = all_fail_bumps.find(team_idx);
        if (fail_bumps_iter == all_fail_bumps.end()) continue;
        fail_bumps = &(fail_bumps_iter->second);
        pt = &vpt[team_idx];

        a_star_reroute_per_team(
            rc,
            *pt,
            *fail_bumps,
            reroute_fail_nets,
            escape_results,
            fail_nets_mutex,
            escape_results_mutex,
            rc_mutex,
            team_idx
        );
    }

    RerouteResult result;
    result.reroute_fail_nets = reroute_fail_nets;
    result.escape_results = escape_results;
    return result;
}
