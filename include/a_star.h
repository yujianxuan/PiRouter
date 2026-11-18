#ifndef A_STAR_H
#define A_STAR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <climits>
#include <algorithm>
#include "json.hpp"
#include "elements.h"
#include "kmean.h"
#include "layer_assignment.h"
#include "encode_coord.h"

using json = nlohmann::json;

struct AStarNode {
    int x, y, layer;
    int g;
    int h;
    int f;
    AStarNode* parent;
    AStarNode(int x_, int y_, int layer_);
    bool operator>(const AStarNode& other) const;
};

struct AStarResult {
    std::vector<AStarNode> path;
    bool hit_node_limit = false;
};

struct EscapeOutput {
    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> escape_results;
    std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>> direct_net_paths;
    std::unordered_map<int, std::vector<std::pair<size_t, int>>> all_fail_bumps;
};

int grid_distance(int x1, int y1, int x2, int y2);

bool is_point_in_convex_polygon(double x, double y, const std::vector<Point>& edges);

bool is_same_direction(int dx1, int dy1, int dl1, int dx2, int dy2, int dl2);

std::vector<std::pair<ResultPoint, ResultPoint>> split_to_segments(const std::vector<AStarNode>& path);

EscapeOutput a_star_escape_with_threads(
    RoutingCase& rc,
    const std::vector<PairTeam>& vpt,
    std::vector<bool> intersect_flags,
    int escape_threads
);

#endif
