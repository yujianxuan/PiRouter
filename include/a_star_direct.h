#ifndef A_STAR_DIRECT_H
#define A_STAR_DIRECT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <climits>
#include <algorithm>
#include "json.hpp"
#include "elements.h"
#include "kmean.h"
#include "layer_assignment.h"
#include "encode_coord.h"

using json = nlohmann::json;

struct AStarNodeDirect {
    int x, y, layer;
    int g;
    int h;
    int f;
    AStarNodeDirect* parent;

    AStarNodeDirect(int x_, int y_, int layer_);

    bool operator>(const AStarNodeDirect& other) const;
};

struct AStarResultDirect {
    std::vector<AStarNodeDirect> path;
};

struct RerouteResult {
    std::vector<std::string> reroute_fail_nets;
    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> escape_results;
};

// 删除 get_coord_node_key_direct 和 get_node_key_direct 声明

int grid_distance_direct(int x1, int y1, int l1, int x2, int y2, int l2);

std::vector<AStarNodeDirect> reconstruct_path_direct(AStarNodeDirect* end_node,
                                                     OccGrid& occupied_set,
                                                     bool commit_path = true);

AStarResultDirect a_star_find_path_direct(const Bump& start,
                                          int max_width, int max_height, int layer_count,
                                          OccGrid& occupied_set,
                                          bool commit_path = true);

// 兜底连接：从任意起点搜索到任意终点（不锁定层），用于 escape+reroute 都失败的 bump
// 把孤立 bump 连到本 net 已布路径的开放端点（partner 的 exit_terminal）
// self_cells：本 net 自己已占的格点，允许穿行（同 net 自重叠合法），解决 bump via stub 把起点困住的问题
AStarResultDirect a_star_connect_points(int sx, int sy, int sl,
                                        int ex, int ey, int el,
                                        int max_width, int max_height, int layer_count,
                                        OccGrid& occupied_set,
                                        const std::unordered_set<uint64_t>& self_cells);

json coord_to_json_direct(int x, int y, int layer);

bool is_same_direction_direct(int dx1, int dy1, int dl1,
                              int dx2, int dy2, int dl2);

std::vector<std::pair<ResultPoint, ResultPoint>> split_to_segments_direct(const std::vector<AStarNodeDirect>& path);

RerouteResult a_star_reroute(
    RoutingCase& rc,
    const std::vector<PairTeam>& vpt,
    std::unordered_map<int, std::vector<std::pair<size_t, int>>> all_fail_bumps
);

RerouteResult a_star_reroute_multithread(
    RoutingCase& rc,
    const std::vector<PairTeam>& vpt,
    std::unordered_map<int, std::vector<std::pair<size_t, int>>> all_fail_bumps
);

#endif
