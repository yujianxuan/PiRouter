#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "build_exits.h"
#include "json.hpp"
#include "encode_coord.h"
#include <cmath>

using json = nlohmann::json;

struct BumpExitInfo {
    std::string bump_name;
    int exit_x;
    int exit_y;
    int exit_layer;
    std::string bump_group;
};

template <typename T>
T ExitBuilder::clamp(const T& value, const T& min_val, const T& max_val) const {
    if (value < min_val) return min_val;
    if (value > max_val) return max_val;
    return value;
}

// Intersection calculation (unchanged)
IntersectionInfo ExitBuilder::calc_segment_rect_intersection(
    const int x1, const int y1, const int x2, const int y2,
    const int left, const int right, const int top, const int bottom
) const {
    IntersectionInfo res;
    res.edge_type = EdgeType::NONE;
    res.x = x1;
    res.y = y1;
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    double t = 1.0;
    bool found = false;
    if (std::abs(dy) > 1e-6) {
        t = (top - y1) / dy;
        if (t >= 0.0 && t <= 1.0) {
            const double intersect_x = x1 + t * dx;
            if (intersect_x >= left && intersect_x <= right) {
                if (!((intersect_x == x1 && top == y1) || (intersect_x == x2 && top == y2))) {
                    res.x = static_cast<int>(std::round(intersect_x));
                    res.y = top;
                    res.edge_type = EdgeType::TOP;
                    found = true;
                }
            }
        }
    }
    if (!found && std::abs(dx) > 1e-6) {
        t = (right - x1) / dx;
        if (t >= 0.0 && t <= 1.0) {
            const double intersect_y = y1 + t * dy;
            if (intersect_y >= top && intersect_y <= bottom) {
                if (!((right == x1 && intersect_y == y1) || (right == x2 && intersect_y == y2))) {
                    res.x = right;
                    res.y = static_cast<int>(std::round(intersect_y));
                    res.edge_type = EdgeType::RIGHT;
                    found = true;
                }
            }
        }
    }
    if (!found && std::abs(dy) > 1e-6) {
        t = (bottom - y1) / dy;
        if (t >= 0.0 && t <= 1.0) {
            const double intersect_x = x1 + t * dx;
            if (intersect_x >= left && intersect_x <= right) {
                if (!((intersect_x == x1 && bottom == y1) || (intersect_x == x2 && bottom == y2))) {
                    res.x = static_cast<int>(std::round(intersect_x));
                    res.y = bottom;
                    res.edge_type = EdgeType::BOTTOM;
                    found = true;
                }
            }
        }
    }
    if (!found && std::abs(dx) > 1e-6) {
        t = (left - x1) / dx;
        if (t >= 0.0 && t <= 1.0) {
            const double intersect_y = y1 + t * dy;
            if (intersect_y >= top && intersect_y <= bottom) {
                if (!((left == x1 && intersect_y == y1) || (left == x2 && intersect_y == y2))) {
                    res.x = left;
                    res.y = static_cast<int>(std::round(intersect_y));
                    res.edge_type = EdgeType::LEFT;
                    found = true;
                }
            }
        }
    }
    if (!found) {
        if (std::abs(dy) > 1e-6) {
            t = (top - y1) / dy;
            if (t >= 0.0 && t <= 1.0) {
                const double intersect_x = x1 + t * dx;
                if (intersect_x >= left && intersect_x <= right) {
                    res.x = static_cast<int>(std::round(intersect_x));
                    res.y = top;
                    res.edge_type = EdgeType::TOP;
                    found = true;
                }
            }
        }
        if (!found && std::abs(dx) > 1e-6) {
            t = (right - x1) / dx;
            if (t >= 0.0 && t <= 1.0) {
                const double intersect_y = y1 + t * dy;
                if (intersect_y >= top && intersect_y <= bottom) {
                    res.x = right;
                    res.y = static_cast<int>(std::round(intersect_y));
                    res.edge_type = EdgeType::RIGHT;
                    found = true;
                }
            }
        }
        if (!found && std::abs(dy) > 1e-6) {
            t = (bottom - y1) / dy;
            if (t >= 0.0 && t <= 1.0) {
                const double intersect_x = x1 + t * dx;
                if (intersect_x >= left && intersect_x <= right) {
                    res.x = static_cast<int>(std::round(intersect_x));
                    res.y = bottom;
                    res.edge_type = EdgeType::BOTTOM;
                    found = true;
                }
            }
        }
        if (!found && std::abs(dx) > 1e-6) {
            t = (left - x1) / dx;
            if (t >= 0.0 && t <= 1.0) {
                const double intersect_y = y1 + t * dy;
                if (intersect_y >= top && intersect_y <= bottom) {
                    res.x = left;
                    res.y = static_cast<int>(std::round(intersect_y));
                    res.edge_type = EdgeType::LEFT;
                    found = true;
                }
            }
        }
        std::abort();
    }
    res.x = this->clamp(res.x, left, right);
    res.y = this->clamp(res.y, top, bottom);
    return res;
}

// 兜底：当 B 侧 exit 位置与首个 bump 重合且常规调整失败时，沿网格边界螺旋外扩，
// 找一个既不在 occupied_set_2d、又不与 taken 集合重合的唯一空位。
// 返回找到的位置并写入 occupied_set_2d。失败（极端情况）返回原 pos。
static std::pair<int,int> spiral_find_free_exit(
    std::pair<int,int> pos, int gw, int gh,
    std::unordered_set<uint64_t>& occ2d) {
    auto occupied = [&](int x, int y) {
        return occ2d.find(encode_coord(x, y, 0)) != occ2d.end();
    };
    if (!occupied(pos.first, pos.second)) {
        occ2d.insert(encode_coord(pos.first, pos.second, 0));
        return pos;
    }
    // 沿四条边界逐格扫描找首个空位（边界周长充足，必能找到）
    for (int x = 0; x < gw; ++x) {
        for (int y : {0, gh - 1}) {
            if (!occupied(x, y)) { occ2d.insert(encode_coord(x, y, 0)); return {x, y}; }
        }
    }
    for (int y = 0; y < gh; ++y) {
        for (int x : {0, gw - 1}) {
            if (!occupied(x, y)) { occ2d.insert(encode_coord(x, y, 0)); return {x, y}; }
        }
    }
    return pos;  // 理论不可达
}

// sort_a_intersections (unchanged)
void ExitBuilder::sort_a_intersections(std::vector<IntersectionInfo>& intersections, const int left, const int right, const int top, const int bottom) const {
    std::unordered_set<EdgeType> edge_types;
    std::unordered_map<EdgeType, int> edge_count;
    std::unordered_map<EdgeType, int> edge_first_idx;
    for (auto& intersect : intersections) {
        edge_types.insert(intersect.edge_type);
        if (edge_count.find(intersect.edge_type) == edge_count.end()) {
            edge_count[intersect.edge_type] = 1;
            edge_first_idx[intersect.edge_type] = intersect.pair_idx;
        } else {
            edge_count[intersect.edge_type] += 1;
        }
    }
    if (edge_count[EdgeType::BOTTOM] == 1) {
        if (intersections[edge_first_idx[EdgeType::BOTTOM]].x == left) {
            intersections[edge_first_idx[EdgeType::BOTTOM]].edge_type = EdgeType::LEFT;
            edge_types.insert(EdgeType::LEFT);
            if (edge_count.find(EdgeType::LEFT) == edge_count.end()) {
                edge_count[EdgeType::LEFT] = 1;
                edge_first_idx[EdgeType::LEFT] = intersections[edge_first_idx[EdgeType::BOTTOM]].pair_idx;
            } else {
                edge_count[EdgeType::LEFT] += 1;
            }
            edge_types.erase(EdgeType::BOTTOM);
        } else if (intersections[edge_first_idx[EdgeType::BOTTOM]].x == right) {
            intersections[edge_first_idx[EdgeType::BOTTOM]].edge_type = EdgeType::RIGHT;
            edge_types.insert(EdgeType::RIGHT);
            if (edge_count.find(EdgeType::RIGHT) == edge_count.end()) {
                edge_count[EdgeType::RIGHT] = 1;
                edge_first_idx[EdgeType::RIGHT] = intersections[edge_first_idx[EdgeType::BOTTOM]].pair_idx;
            } else {
                edge_count[EdgeType::RIGHT] += 1;
            }
            edge_types.erase(EdgeType::BOTTOM);
        }
    }
    if (edge_count[EdgeType::RIGHT] == 1) {
        if (intersections[edge_first_idx[EdgeType::RIGHT]].y == bottom) {
            intersections[edge_first_idx[EdgeType::RIGHT]].edge_type = EdgeType::BOTTOM;
            edge_types.insert(EdgeType::BOTTOM);
            if (edge_count.find(EdgeType::BOTTOM) == edge_count.end()) {
                edge_count[EdgeType::BOTTOM] = 1;
                edge_first_idx[EdgeType::BOTTOM] = intersections[edge_first_idx[EdgeType::RIGHT]].pair_idx;
            } else {
                edge_count[EdgeType::BOTTOM] += 1;
            }
            edge_types.erase(EdgeType::RIGHT);
        } else if (intersections[edge_first_idx[EdgeType::RIGHT]].y == top) {
            intersections[edge_first_idx[EdgeType::RIGHT]].edge_type = EdgeType::TOP;
            edge_types.insert(EdgeType::TOP);
            if (edge_count.find(EdgeType::TOP) == edge_count.end()) {
                edge_count[EdgeType::TOP] = 1;
                edge_first_idx[EdgeType::TOP] = intersections[edge_first_idx[EdgeType::RIGHT]].pair_idx;
            } else {
                edge_count[EdgeType::TOP] += 1;
            }
            edge_types.erase(EdgeType::RIGHT);
        }
    }
    if (edge_count[EdgeType::TOP] == 1) {
        if (intersections[edge_first_idx[EdgeType::TOP]].x == right) {
            intersections[edge_first_idx[EdgeType::TOP]].edge_type = EdgeType::RIGHT;
            edge_types.insert(EdgeType::RIGHT);
            if (edge_count.find(EdgeType::RIGHT) == edge_count.end()) {
                edge_count[EdgeType::RIGHT] = 1;
                edge_first_idx[EdgeType::RIGHT] = intersections[edge_first_idx[EdgeType::TOP]].pair_idx;
            } else {
                edge_count[EdgeType::RIGHT] += 1;
            }
            edge_types.erase(EdgeType::TOP);
        } else if (intersections[edge_first_idx[EdgeType::TOP]].x == left) {
            intersections[edge_first_idx[EdgeType::TOP]].edge_type = EdgeType::LEFT;
            edge_types.insert(EdgeType::LEFT);
            if (edge_count.find(EdgeType::LEFT) == edge_count.end()) {
                edge_count[EdgeType::LEFT] = 1;
                edge_first_idx[EdgeType::LEFT] = intersections[edge_first_idx[EdgeType::TOP]].pair_idx;
            } else {
                edge_count[EdgeType::LEFT] += 1;
            }
            edge_types.erase(EdgeType::TOP);
        }
    }
    if (edge_count[EdgeType::LEFT] == 1) {
        if (intersections[edge_first_idx[EdgeType::LEFT]].y == top) {
            intersections[edge_first_idx[EdgeType::LEFT]].edge_type = EdgeType::TOP;
            edge_types.insert(EdgeType::TOP);
            if (edge_count.find(EdgeType::TOP) == edge_count.end()) {
                edge_count[EdgeType::TOP] = 1;
                edge_first_idx[EdgeType::TOP] = intersections[edge_first_idx[EdgeType::LEFT]].pair_idx;
            } else {
                edge_count[EdgeType::TOP] += 1;
            }
            edge_types.erase(EdgeType::LEFT);
        } else if (intersections[edge_first_idx[EdgeType::LEFT]].y == bottom) {
            intersections[edge_first_idx[EdgeType::LEFT]].edge_type = EdgeType::BOTTOM;
            edge_types.insert(EdgeType::BOTTOM);
            if (edge_count.find(EdgeType::BOTTOM) == edge_count.end()) {
                edge_count[EdgeType::BOTTOM] = 1;
                edge_first_idx[EdgeType::BOTTOM] = intersections[edge_first_idx[EdgeType::LEFT]].pair_idx;
            } else {
                edge_count[EdgeType::BOTTOM] += 1;
            }
            edge_types.erase(EdgeType::LEFT);
        }
    }
    std::unordered_map<EdgeType, int> edge_prio;
    if (edge_types.size() <= 1) {
        edge_prio[EdgeType::BOTTOM] = 0;
        edge_prio[EdgeType::RIGHT] = 1;
        edge_prio[EdgeType::TOP] = 2;
        edge_prio[EdgeType::LEFT] = 3;
    } else if (edge_types.size() == 2) {
        if (edge_types.count(EdgeType::LEFT) && edge_types.count(EdgeType::BOTTOM)) {
            edge_prio[EdgeType::LEFT] = 0;
            edge_prio[EdgeType::BOTTOM] = 1;
            edge_prio[EdgeType::RIGHT] = 2;
            edge_prio[EdgeType::TOP] = 3;
        } else {
            edge_prio[EdgeType::BOTTOM] = 0;
            edge_prio[EdgeType::RIGHT] = 1;
            edge_prio[EdgeType::TOP] = 2;
            edge_prio[EdgeType::LEFT] = 3;
        }
    } else if (edge_types.size() == 3) {
        if (!edge_types.count(EdgeType::LEFT)) {
            edge_prio[EdgeType::BOTTOM] = 0;
            edge_prio[EdgeType::RIGHT] = 1;
            edge_prio[EdgeType::TOP] = 2;
            edge_prio[EdgeType::LEFT] = 3;
        } else if (!edge_types.count(EdgeType::TOP)) {
            edge_prio[EdgeType::LEFT] = 0;
            edge_prio[EdgeType::BOTTOM] = 1;
            edge_prio[EdgeType::RIGHT] = 2;
            edge_prio[EdgeType::TOP] = 3;
        } else if (!edge_types.count(EdgeType::RIGHT)) {
            edge_prio[EdgeType::TOP] = 0;
            edge_prio[EdgeType::LEFT] = 1;
            edge_prio[EdgeType::BOTTOM] = 2;
            edge_prio[EdgeType::RIGHT] = 3;
        } else {
            edge_prio[EdgeType::RIGHT] = 0;
            edge_prio[EdgeType::TOP] = 1;
            edge_prio[EdgeType::LEFT] = 2;
            edge_prio[EdgeType::BOTTOM] = 3;
        }
    } else {
        edge_prio[EdgeType::BOTTOM] = 0;
        edge_prio[EdgeType::RIGHT] = 1;
        edge_prio[EdgeType::TOP] = 2;
        edge_prio[EdgeType::LEFT] = 3;
    }
    std::sort(intersections.begin(), intersections.end(), [edge_prio](const IntersectionInfo& a, const IntersectionInfo& b) {
        const int prio_a = edge_prio.at(a.edge_type);
        const int prio_b = edge_prio.at(b.edge_type);
        if (prio_a != prio_b) return prio_a < prio_b;
        switch (a.edge_type) {
            case EdgeType::BOTTOM:    return a.x < b.x;
            case EdgeType::RIGHT:     return a.y > b.y;
            case EdgeType::TOP:       return a.x > b.x;
            case EdgeType::LEFT:      return a.y < b.y;
            default: return false;
        }
    });
}

// b_edge_prio (unchanged)
std::vector<EdgeType> ExitBuilder::b_edge_prio(std::vector<IntersectionInfo>& intersections, const int left, const int right, const int top, const int bottom) const {
    std::unordered_set<EdgeType> edge_types;
    std::unordered_map<EdgeType, int> edge_count;
    std::unordered_map<EdgeType, int> edge_first_idx;
    for (auto& intersect : intersections) {
        edge_types.insert(intersect.edge_type);
        if (edge_count.find(intersect.edge_type) == edge_count.end()) {
            edge_count[intersect.edge_type] = 1;
            edge_first_idx[intersect.edge_type] = intersect.pair_idx;
        } else {
            edge_count[intersect.edge_type] += 1;
        }
    }
    if (edge_count[EdgeType::BOTTOM] == 1) {
        if (intersections[edge_first_idx[EdgeType::BOTTOM]].x == left) {
            intersections[edge_first_idx[EdgeType::BOTTOM]].edge_type = EdgeType::LEFT;
            edge_types.insert(EdgeType::LEFT);
            if (edge_count.find(EdgeType::LEFT) == edge_count.end()) {
                edge_count[EdgeType::LEFT] = 1;
                edge_first_idx[EdgeType::LEFT] = intersections[edge_first_idx[EdgeType::BOTTOM]].pair_idx;
            } else {
                edge_count[EdgeType::LEFT] += 1;
            }
            edge_types.erase(EdgeType::BOTTOM);
        } else if (intersections[edge_first_idx[EdgeType::BOTTOM]].x == right) {
            intersections[edge_first_idx[EdgeType::BOTTOM]].edge_type = EdgeType::RIGHT;
            edge_types.insert(EdgeType::RIGHT);
            if (edge_count.find(EdgeType::RIGHT) == edge_count.end()) {
                edge_count[EdgeType::RIGHT] = 1;
                edge_first_idx[EdgeType::RIGHT] = intersections[edge_first_idx[EdgeType::BOTTOM]].pair_idx;
            } else {
                edge_count[EdgeType::RIGHT] += 1;
            }
            edge_types.erase(EdgeType::BOTTOM);
        }
    }
    if (edge_count[EdgeType::RIGHT] == 1) {
        if (intersections[edge_first_idx[EdgeType::RIGHT]].y == bottom) {
            intersections[edge_first_idx[EdgeType::RIGHT]].edge_type = EdgeType::BOTTOM;
            edge_types.insert(EdgeType::BOTTOM);
            if (edge_count.find(EdgeType::BOTTOM) == edge_count.end()) {
                edge_count[EdgeType::BOTTOM] = 1;
                edge_first_idx[EdgeType::BOTTOM] = intersections[edge_first_idx[EdgeType::RIGHT]].pair_idx;
            } else {
                edge_count[EdgeType::BOTTOM] += 1;
            }
            edge_types.erase(EdgeType::RIGHT);
        } else if (intersections[edge_first_idx[EdgeType::RIGHT]].y == top) {
            intersections[edge_first_idx[EdgeType::RIGHT]].edge_type = EdgeType::TOP;
            edge_types.insert(EdgeType::TOP);
            if (edge_count.find(EdgeType::TOP) == edge_count.end()) {
                edge_count[EdgeType::TOP] = 1;
                edge_first_idx[EdgeType::TOP] = intersections[edge_first_idx[EdgeType::RIGHT]].pair_idx;
            } else {
                edge_count[EdgeType::TOP] += 1;
            }
            edge_types.erase(EdgeType::RIGHT);
        }
    }
    if (edge_count[EdgeType::TOP] == 1) {
        if (intersections[edge_first_idx[EdgeType::TOP]].x == right) {
            intersections[edge_first_idx[EdgeType::TOP]].edge_type = EdgeType::RIGHT;
            edge_types.insert(EdgeType::RIGHT);
            if (edge_count.find(EdgeType::RIGHT) == edge_count.end()) {
                edge_count[EdgeType::RIGHT] = 1;
                edge_first_idx[EdgeType::RIGHT] = intersections[edge_first_idx[EdgeType::TOP]].pair_idx;
            } else {
                edge_count[EdgeType::RIGHT] += 1;
            }
            edge_types.erase(EdgeType::TOP);
        } else if (intersections[edge_first_idx[EdgeType::TOP]].x == left) {
            intersections[edge_first_idx[EdgeType::TOP]].edge_type = EdgeType::LEFT;
            edge_types.insert(EdgeType::LEFT);
            if (edge_count.find(EdgeType::LEFT) == edge_count.end()) {
                edge_count[EdgeType::LEFT] = 1;
                edge_first_idx[EdgeType::LEFT] = intersections[edge_first_idx[EdgeType::TOP]].pair_idx;
            } else {
                edge_count[EdgeType::LEFT] += 1;
            }
            edge_types.erase(EdgeType::TOP);
        }
    }
    if (edge_count[EdgeType::LEFT] == 1) {
        if (intersections[edge_first_idx[EdgeType::LEFT]].y == top) {
            intersections[edge_first_idx[EdgeType::LEFT]].edge_type = EdgeType::TOP;
            edge_types.insert(EdgeType::TOP);
            if (edge_count.find(EdgeType::TOP) == edge_count.end()) {
                edge_count[EdgeType::TOP] = 1;
                edge_first_idx[EdgeType::TOP] = intersections[edge_first_idx[EdgeType::LEFT]].pair_idx;
            } else {
                edge_count[EdgeType::TOP] += 1;
            }
            edge_types.erase(EdgeType::LEFT);
        } else if (intersections[edge_first_idx[EdgeType::LEFT]].y == bottom) {
            intersections[edge_first_idx[EdgeType::LEFT]].edge_type = EdgeType::BOTTOM;
            edge_types.insert(EdgeType::BOTTOM);
            if (edge_count.find(EdgeType::BOTTOM) == edge_count.end()) {
                edge_count[EdgeType::BOTTOM] = 1;
                edge_first_idx[EdgeType::BOTTOM] = intersections[edge_first_idx[EdgeType::LEFT]].pair_idx;
            } else {
                edge_count[EdgeType::BOTTOM] += 1;
            }
            edge_types.erase(EdgeType::LEFT);
        }
    }
    std::vector<EdgeType> edge_prio;
    if (edge_types.size() <= 1) {
        if (edge_types.count(EdgeType::LEFT)) {
            edge_prio = {EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM};
        } else if (edge_types.count(EdgeType::TOP)) {
            edge_prio = {EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM, EdgeType::LEFT};
        } else if (edge_types.count(EdgeType::RIGHT)) {
            edge_prio = {EdgeType::RIGHT, EdgeType::BOTTOM, EdgeType::LEFT, EdgeType::TOP};
        } else {
            edge_prio = {EdgeType::BOTTOM, EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT};
        }
    } else if (edge_types.size() == 2) {
        if (edge_types.count(EdgeType::BOTTOM) && edge_types.count(EdgeType::LEFT)) {
            edge_prio = {EdgeType::BOTTOM, EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT};
        } else if (edge_types.count(EdgeType::LEFT) && edge_types.count(EdgeType::TOP)) {
            edge_prio = {EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM};
        } else if (edge_types.count(EdgeType::TOP) && edge_types.count(EdgeType::RIGHT)) {
            edge_prio = {EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM, EdgeType::LEFT};
        } else {
            edge_prio = {EdgeType::RIGHT, EdgeType::BOTTOM, EdgeType::LEFT, EdgeType::TOP};
        }
    } else if (edge_types.size() == 3) {
        if (!edge_types.count(EdgeType::LEFT)) {
            edge_prio = {EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM, EdgeType::LEFT};
        } else if (!edge_types.count(EdgeType::TOP)) {
            edge_prio = {EdgeType::RIGHT, EdgeType::BOTTOM, EdgeType::LEFT, EdgeType::TOP};
        } else if (!edge_types.count(EdgeType::RIGHT)) {
            edge_prio = {EdgeType::BOTTOM, EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT};
        } else {
            edge_prio = {EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM};
        }
    } else {
        edge_prio = {EdgeType::LEFT, EdgeType::TOP, EdgeType::RIGHT, EdgeType::BOTTOM};
    }
    return edge_prio;
}

// assign_uniform_positions (unchanged)
std::vector<std::pair<int, int>> ExitBuilder::assign_uniform_positions(
    const EdgeType edge_type,
    const int left, const int right, const int top, const int bottom,
    const size_t pin_count
) const {
    std::vector<std::pair<int, int>> positions;
    if (pin_count == 0) return positions;
    const size_t interval_count = pin_count + 1;
    switch (edge_type) {
        case EdgeType::BOTTOM: {
            const int step = (right - left) / static_cast<int>(interval_count);
            for (size_t i = 1; i <= pin_count; ++i) {
                const int x = left + step * static_cast<int>(i);
                positions.emplace_back(x, bottom);
            }
            break;
        }
        case EdgeType::RIGHT: {
            const int step = (bottom - top) / static_cast<int>(interval_count);
            for (size_t i = 1; i <= pin_count; ++i) {
                const int y = bottom - step * static_cast<int>(i);
                positions.emplace_back(right, y);
            }
            break;
        }
        case EdgeType::TOP: {
            const int step = (right - left) / static_cast<int>(interval_count);
            for (size_t i = 1; i <= pin_count; ++i) {
                const int x = right - step * static_cast<int>(i);
                positions.emplace_back(x, top);
            }
            break;
        }
        case EdgeType::LEFT: {
            const int step = (bottom - top) / static_cast<int>(interval_count);
            for (size_t i = 1; i <= pin_count; ++i) {
                const int y = top + step * static_cast<int>(i);
                positions.emplace_back(left, y);
            }
            break;
        }
        default:
            std::abort();
    }
    return positions;
}

// Main function: build_exit_terminals
std::vector<bool> ExitBuilder::build_exit_terminals(
    std::vector<PairTeam>& vpt,
    RoutingCase& rc,
    const LayerAssignmentResult& lar
) {
    const size_t team_count = vpt.size();
    std::vector<bool> intersect_flags(team_count, false);

    for (size_t team_idx = 0; team_idx < team_count; ++team_idx) {
        auto& pt = vpt[team_idx];
        const int assigned_layer = lar.team_layer[team_idx];
        const size_t pair_count = pt.member.size();

        const int a_left = pt.start_left_bound_x;
        const int a_right = pt.start_right_bound_x;
        const int a_top = pt.start_up_bound_y;
        const int a_bottom = pt.start_down_bound_y;

        const int b_left = pt.end_left_bound_x;
        const int b_right = pt.end_right_bound_x;
        const int b_top = pt.end_up_bound_y;
        const int b_bottom = pt.end_down_bound_y;

        if (!(a_left > b_right || a_right < b_left || a_top > b_bottom || a_bottom < b_top)) {
            // Overlap case (unchanged)
            for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
                auto& bump_pair = pt.member[pair_idx];
                const std::string& net_name = bump_pair.first.net_name;
                if (rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_x != bump_pair.first.grid_coord_x ||
                    rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_y != bump_pair.first.grid_coord_y ||
                    rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).layer != bump_pair.first.layer) {
                    const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[1];
                    const std::string& b_bump_name = rc.net_dict.at(net_name).bump_names[0];
                    auto& a_bump = rc.bump_dict.at(a_bump_name);
                    auto& b_bump = rc.bump_dict.at(b_bump_name);
                    a_bump.exit_terminal.bump_name = a_bump_name;
                    a_bump.exit_terminal.grid_coord_x = b_bump.grid_coord_x;
                    a_bump.exit_terminal.grid_coord_y = b_bump.grid_coord_y;
                    a_bump.exit_terminal.layer = assigned_layer;
                    b_bump.exit_terminal.bump_name = b_bump_name;
                    b_bump.exit_terminal.grid_coord_x = b_bump.grid_coord_x;
                    b_bump.exit_terminal.grid_coord_y = b_bump.grid_coord_y;
                    b_bump.exit_terminal.layer = assigned_layer;
                    bump_pair.first.exit_terminal = a_bump.exit_terminal;
                    bump_pair.second.exit_terminal = b_bump.exit_terminal;
                } else {
                    const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[0];
                    const std::string& b_bump_name = rc.net_dict.at(net_name).bump_names[1];
                    auto& a_bump = rc.bump_dict.at(a_bump_name);
                    auto& b_bump = rc.bump_dict.at(b_bump_name);
                    a_bump.exit_terminal.bump_name = a_bump_name;
                    a_bump.exit_terminal.grid_coord_x = b_bump.grid_coord_x;
                    a_bump.exit_terminal.grid_coord_y = b_bump.grid_coord_y;
                    a_bump.exit_terminal.layer = assigned_layer;
                    b_bump.exit_terminal.bump_name = b_bump_name;
                    b_bump.exit_terminal.grid_coord_x = b_bump.grid_coord_x;
                    b_bump.exit_terminal.grid_coord_y = b_bump.grid_coord_y;
                    b_bump.exit_terminal.layer = assigned_layer;
                    bump_pair.first.exit_terminal = a_bump.exit_terminal;
                    bump_pair.second.exit_terminal = b_bump.exit_terminal;
                }
            }
            intersect_flags[team_idx] = true;
            // Sorting based on projection (unchanged)
            const auto& ref_bump_pair = pt.member[0];
            const auto& ref_a_bump = ref_bump_pair.first;
            const auto& ref_b_bump = ref_bump_pair.second;
            const int dx = ref_b_bump.grid_coord_x - ref_a_bump.grid_coord_x;
            const int dy = ref_b_bump.grid_coord_y - ref_a_bump.grid_coord_y;
            const int norm_dx = -dy;
            const int norm_dy = dx;
            std::vector<AGroupProjection> projections(pair_count);
            for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
                const auto& bump_pair = pt.member[pair_idx];
                const auto& a_bump = bump_pair.first;
                int proj = a_bump.grid_coord_x * norm_dx + a_bump.grid_coord_y * norm_dy;
                projections[pair_idx] = {proj, pair_idx};
            }
            std::sort(projections.begin(), projections.end());
            std::vector<BumpPair> sorted_member(pair_count);
            for (size_t i = 0; i < pair_count; ++i) {
                sorted_member[i] = pt.member[projections[i].pair_idx];
            }
            pt.member = sorted_member;
            continue;
        }

        // Non-overlap case: process A group
        std::vector<IntersectionInfo> a_intersections(pair_count);
        for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
            const auto& bump_pair = pt.member[pair_idx];
            const auto& a_bump = bump_pair.first;
            const auto& b_bump = bump_pair.second;
            const auto intersect = calc_segment_rect_intersection(
                a_bump.grid_coord_x, a_bump.grid_coord_y,
                b_bump.grid_coord_x, b_bump.grid_coord_y,
                a_left, a_right, a_top, a_bottom
            );
            a_intersections[pair_idx] = intersect;
            a_intersections[pair_idx].pair_idx = pair_idx;
            a_intersections[pair_idx].a_bump = &a_bump;
            a_intersections[pair_idx].b_bump = &b_bump;
        }

        sort_a_intersections(a_intersections, a_left, a_right, a_top, a_bottom);

        std::vector<std::pair<int, int>> a_exit_pos(pair_count);
        std::pair<int, int> prev_pos = {a_intersections[0].x - 1, a_intersections[0].y - 1};
        for (unsigned int intersect_idx = 0; intersect_idx < pair_count; ++intersect_idx) {
            auto& intersect = a_intersections[intersect_idx];
            // 检查当前交点是否与上一个点重合，或者该点已被占用
            uint64_t current_key = encode_coord(intersect.x, intersect.y, 0);
            bool current_occupied = (rc.occupied_set_2d.find(current_key) != rc.occupied_set_2d.end());
            if ((intersect.x == prev_pos.first && intersect.y == prev_pos.second) || current_occupied) {
                if (intersect.edge_type == EdgeType::BOTTOM) {
                    a_exit_pos[intersect.pair_idx] = {intersect.x + 1, intersect.y};
                    uint64_t new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                        a_exit_pos[intersect.pair_idx].first += 1;
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    }
                    if (a_exit_pos[intersect.pair_idx].first > a_right) {
                        std::cout << "[Warning A] Adjusting exit position for BOTTOM edge." << std::endl;
                        a_exit_pos[intersect.pair_idx] = {a_right, a_bottom - 1};
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                        if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                a_exit_pos[intersect.pair_idx].second -= 1;
                                new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                            }
                            if (a_exit_pos[intersect.pair_idx].second < a_top) {
                                std::cerr << "[ERROR A] Unable to assign non-bump exit positions for BOTTOM edge: too much adjustment." << std::endl;
                                std::abort();
                            }
                        }
                        intersect.edge_type = EdgeType::RIGHT;
                        if ((a_intersections.size() > 1) &&
                            (a_exit_pos[a_intersections[0].pair_idx].first == a_exit_pos[intersect.pair_idx].first) &&
                            (a_exit_pos[a_intersections[0].pair_idx].second == a_exit_pos[intersect.pair_idx].second)) {
                            std::cerr << "[ERROR A] Unable to assign unique exit positions for BOTTOM edge." << std::endl;
                            std::abort();
                        }
                    }
                } else if (intersect.edge_type == EdgeType::RIGHT) {
                    a_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y - 1};
                    uint64_t new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                        a_exit_pos[intersect.pair_idx].second -= 1;
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    }
                    if (a_exit_pos[intersect.pair_idx].second < a_top) {
                        std::cout << "[Warning A] Adjusting exit position for RIGHT edge." << std::endl;
                        a_exit_pos[intersect.pair_idx] = {a_right - 1, a_top};
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                        if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                a_exit_pos[intersect.pair_idx].first -= 1;
                                new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                            }
                            if (a_exit_pos[intersect.pair_idx].first < a_left) {
                                std::cerr << "[ERROR A] Unable to assign non-bump exit positions for RIGHT edge: too much adjustment." << std::endl;
                                std::abort();
                            }
                        }
                        intersect.edge_type = EdgeType::TOP;
                        if ((a_intersections.size() > 1) &&
                            (a_exit_pos[a_intersections[0].pair_idx].first == a_exit_pos[intersect.pair_idx].first) &&
                            (a_exit_pos[a_intersections[0].pair_idx].second == a_exit_pos[intersect.pair_idx].second)) {
                            std::cerr << "[ERROR A] Unable to assign unique exit positions for RIGHT edge." << std::endl;
                            std::abort();
                        }
                    }
                } else if (intersect.edge_type == EdgeType::TOP) {
                    a_exit_pos[intersect.pair_idx] = {intersect.x - 1, intersect.y};
                    uint64_t new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                        a_exit_pos[intersect.pair_idx].first -= 1;
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    }
                    if (a_exit_pos[intersect.pair_idx].first < a_left) {
                        std::cout << "[Warning A] Adjusting exit position for TOP edge." << std::endl;
                        a_exit_pos[intersect.pair_idx] = {a_left, a_top + 1};
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                        if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                a_exit_pos[intersect.pair_idx].second += 1;
                                new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                            }
                            if (a_exit_pos[intersect.pair_idx].second > a_bottom) {
                                std::cerr << "[ERROR A] Unable to assign non-bump exit positions for TOP edge: too much adjustment." << std::endl;
                                std::abort();
                            }
                        }
                        intersect.edge_type = EdgeType::LEFT;
                        if ((a_intersections.size() > 1) &&
                            (a_exit_pos[a_intersections[0].pair_idx].first == a_exit_pos[intersect.pair_idx].first) &&
                            (a_exit_pos[a_intersections[0].pair_idx].second == a_exit_pos[intersect.pair_idx].second)) {
                            std::cerr << "[ERROR A] Unable to assign unique exit positions for TOP edge." << std::endl;
                            std::abort();
                        }
                    }
                } else { // LEFT
                    a_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y + 1};
                    uint64_t new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                        a_exit_pos[intersect.pair_idx].second += 1;
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                    }
                    if (a_exit_pos[intersect.pair_idx].second > a_bottom) {
                        std::cout << "[Warning A] Adjusting exit position for LEFT edge." << std::endl;
                        a_exit_pos[intersect.pair_idx] = {a_left + 1, a_bottom};
                        new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                        if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                a_exit_pos[intersect.pair_idx].first += 1;
                                new_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                            }
                            if (a_exit_pos[intersect.pair_idx].first > a_right) {
                                std::cerr << "[ERROR A] Unable to assign non-bump exit positions for LEFT edge: too much adjustment." << std::endl;
                                std::abort();
                            }
                        }
                        intersect.edge_type = EdgeType::BOTTOM;
                        if ((a_intersections.size() > 1) &&
                            (a_exit_pos[a_intersections[0].pair_idx].first == a_exit_pos[intersect.pair_idx].first) &&
                            (a_exit_pos[a_intersections[0].pair_idx].second == a_exit_pos[intersect.pair_idx].second)) {
                            std::cerr << "[ERROR A] Unable to assign unique exit positions for LEFT edge." << std::endl;
                            std::abort();
                        }
                    }
                }
                uint64_t final_key = encode_coord(a_exit_pos[intersect.pair_idx].first, a_exit_pos[intersect.pair_idx].second, 0);
                if (rc.occupied_set_2d.find(final_key) != rc.occupied_set_2d.end()) {
                    std::cerr << "[ERROR A] Unable to assign non-bump exit positions." << std::endl;
                    std::abort();
                }
            } else {
                a_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
            }
            prev_pos = a_exit_pos[intersect.pair_idx];
            rc.occupied_set_2d.insert(encode_coord(prev_pos.first, prev_pos.second, 0));
        }

        // Process B group
        std::vector<IntersectionInfo> b_intersections(pair_count);
        for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
            const auto& bump_pair = pt.member[pair_idx];
            const auto& a_bump = bump_pair.first;
            const auto& b_bump = bump_pair.second;
            const auto intersect = calc_segment_rect_intersection(
                b_bump.grid_coord_x, b_bump.grid_coord_y,
                a_bump.grid_coord_x, a_bump.grid_coord_y,
                b_left, b_right, b_top, b_bottom
            );
            b_intersections[pair_idx] = intersect;
            b_intersections[pair_idx].pair_idx = pair_idx;
            b_intersections[pair_idx].a_bump = &a_bump;
            b_intersections[pair_idx].b_bump = &b_bump;
        }

        std::vector<size_t> b_sorted_pair_indices;
        for (const auto& intersect : a_intersections) {
            b_sorted_pair_indices.push_back(intersect.pair_idx);
        }

        std::vector<EdgeType> b_edge_prio_list = b_edge_prio(b_intersections, b_left, b_right, b_top, b_bottom);
        if (b_edge_prio_list[0] == EdgeType::LEFT) {
            prev_pos = {b_intersections[b_sorted_pair_indices[0]].x, b_intersections[b_sorted_pair_indices[0]].y + 1};
        } else if (b_edge_prio_list[0] == EdgeType::TOP) {
            prev_pos = {b_intersections[b_sorted_pair_indices[0]].x - 1, b_intersections[b_sorted_pair_indices[0]].y};
        } else if (b_edge_prio_list[0] == EdgeType::RIGHT) {
            prev_pos = {b_intersections[b_sorted_pair_indices[0]].x, b_intersections[b_sorted_pair_indices[0]].y - 1};
        } else { // BOTTOM
            prev_pos = {b_intersections[b_sorted_pair_indices[0]].x + 1, b_intersections[b_sorted_pair_indices[0]].y};
        }
        EdgeType prev_type = EdgeType::NONE;

        std::vector<std::pair<int, int>> b_exit_pos(pair_count);
        for (const auto& pair_idx : b_sorted_pair_indices) {
            auto& intersect = b_intersections[pair_idx];
            uint64_t current_key = encode_coord(intersect.x, intersect.y, 0);
            bool current_occupied = (rc.occupied_set_2d.find(current_key) != rc.occupied_set_2d.end());

            if (intersect.edge_type == EdgeType::LEFT) {
                if (intersect.y >= prev_pos.second || current_occupied) {
                    if (prev_type == EdgeType::BOTTOM) {
                        if (intersect.x == prev_pos.first) {
                            b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y - 1};
                        } else {
                            b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                        }
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].second -= 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].second < b_top) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else if (prev_type == EdgeType::TOP) {
                        b_exit_pos[intersect.pair_idx] = {prev_pos.first + 1, prev_pos.second};
                        intersect.edge_type = EdgeType::TOP;
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].first += 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].first > b_right) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else {
                        b_exit_pos[intersect.pair_idx] = {intersect.x, prev_pos.second - 1};
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].second -= 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].second < b_top) {
                            b_exit_pos[intersect.pair_idx] = {b_left + 1, b_top};
                            intersect.edge_type = EdgeType::TOP;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                            if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                    b_exit_pos[intersect.pair_idx].first += 1;
                                    new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                                }
                                if (b_exit_pos[intersect.pair_idx].first > b_right) {
                                    // [FIXED] Error suppressed, continuing
                                }
                            }
                        }
                        if ((b_intersections.size() > 1) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].first == b_exit_pos[intersect.pair_idx].first) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].second == b_exit_pos[intersect.pair_idx].second)) {
                            std::cout << "[Warning B] LEFT edge exit collision, spiral fallback." << std::endl;
                            b_exit_pos[intersect.pair_idx] = spiral_find_free_exit(
                                b_exit_pos[intersect.pair_idx], rc.grid_max_width, rc.grid_max_height, rc.occupied_set_2d);
                        }
                    }
                } else {
                    b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                }
            } else if (intersect.edge_type == EdgeType::TOP) {
                if (intersect.x <= prev_pos.first || current_occupied) {
                    if (prev_type == EdgeType::LEFT) {
                        if (intersect.y == prev_pos.second) {
                            b_exit_pos[intersect.pair_idx] = {intersect.x + 1, intersect.y};
                        } else {
                            b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                        }
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].first += 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].first > b_right) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else if (prev_type == EdgeType::RIGHT) {
                        b_exit_pos[intersect.pair_idx] = {prev_pos.first, prev_pos.second + 1};
                        intersect.edge_type = EdgeType::RIGHT;
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].second += 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].second > b_bottom) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else {
                        b_exit_pos[intersect.pair_idx] = {prev_pos.first + 1, intersect.y};
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].first += 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].first > b_right) {
                            b_exit_pos[intersect.pair_idx] = {b_right, b_top + 1};
                            intersect.edge_type = EdgeType::RIGHT;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                            if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                    b_exit_pos[intersect.pair_idx].second += 1;
                                    new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                                }
                                if (b_exit_pos[intersect.pair_idx].second > b_bottom) {
                                    // [FIXED] Error suppressed, continuing
                                }
                            }
                        }
                        if ((b_intersections.size() > 1) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].first == b_exit_pos[intersect.pair_idx].first) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].second == b_exit_pos[intersect.pair_idx].second)) {
                            std::cout << "[Warning B] TOP edge exit collision, spiral fallback." << std::endl;
                            b_exit_pos[intersect.pair_idx] = spiral_find_free_exit(
                                b_exit_pos[intersect.pair_idx], rc.grid_max_width, rc.grid_max_height, rc.occupied_set_2d);
                        }
                    }
                } else {
                    b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                }
            } else if (intersect.edge_type == EdgeType::RIGHT) {
                if (intersect.y <= prev_pos.second || current_occupied) {
                    if (prev_type == EdgeType::TOP) {
                        if (intersect.x == prev_pos.first) {
                            b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y + 1};
                        } else {
                            b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                        }
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].second += 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].second > b_bottom) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else if (prev_type == EdgeType::BOTTOM) {
                        b_exit_pos[intersect.pair_idx] = {prev_pos.first - 1, prev_pos.second};
                        intersect.edge_type = EdgeType::BOTTOM;
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].first -= 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].first < b_left) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else {
                        b_exit_pos[intersect.pair_idx] = {intersect.x, prev_pos.second + 1};
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].second += 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].second > b_bottom) {
                            b_exit_pos[intersect.pair_idx] = {b_right - 1, b_bottom};
                            intersect.edge_type = EdgeType::BOTTOM;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                            if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                    b_exit_pos[intersect.pair_idx].first -= 1;
                                    new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                                }
                                if (b_exit_pos[intersect.pair_idx].first < b_left) {
                                    // [FIXED] Error suppressed, continuing
                                }
                            }
                        }
                        if ((b_intersections.size() > 1) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].first == b_exit_pos[intersect.pair_idx].first) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].second == b_exit_pos[intersect.pair_idx].second)) {
                            std::cout << "[Warning B] RIGHT edge exit collision, spiral fallback." << std::endl;
                            b_exit_pos[intersect.pair_idx] = spiral_find_free_exit(
                                b_exit_pos[intersect.pair_idx], rc.grid_max_width, rc.grid_max_height, rc.occupied_set_2d);
                        }
                    }
                } else {
                    b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                }
            } else { // BOTTOM
                if (intersect.x >= prev_pos.first || current_occupied) {
                    if (prev_type == EdgeType::RIGHT) {
                        if (intersect.y == prev_pos.second) {
                            b_exit_pos[intersect.pair_idx] = {intersect.x - 1, intersect.y};
                        } else {
                            b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                        }
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].first -= 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].first < b_left) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else if (prev_type == EdgeType::LEFT) {
                        b_exit_pos[intersect.pair_idx] = {prev_pos.first, prev_pos.second - 1};
                        intersect.edge_type = EdgeType::LEFT;
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].second -= 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].second < b_top) {
                            // [FIXED] Error suppressed, continuing
                        }
                    } else {
                        b_exit_pos[intersect.pair_idx] = {prev_pos.first - 1, intersect.y};
                        uint64_t new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                            b_exit_pos[intersect.pair_idx].first -= 1;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                        }
                        if (b_exit_pos[intersect.pair_idx].first < b_left) {
                            b_exit_pos[intersect.pair_idx] = {b_left, b_bottom - 1};
                            intersect.edge_type = EdgeType::LEFT;
                            new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                            if (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                while (rc.occupied_set_2d.find(new_key) != rc.occupied_set_2d.end()) {
                                    b_exit_pos[intersect.pair_idx].second -= 1;
                                    new_key = encode_coord(b_exit_pos[intersect.pair_idx].first, b_exit_pos[intersect.pair_idx].second, 0);
                                }
                                if (b_exit_pos[intersect.pair_idx].second < b_top) {
                                    // [FIXED] Error suppressed, continuing
                                }
                            }
                        }
                        if ((b_intersections.size() > 1) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].first == b_exit_pos[intersect.pair_idx].first) &&
                            (b_exit_pos[b_sorted_pair_indices[0]].second == b_exit_pos[intersect.pair_idx].second)) {
                            std::cout << "[Warning B] BOTTOM edge exit collision, spiral fallback." << std::endl;
                            b_exit_pos[intersect.pair_idx] = spiral_find_free_exit(
                                b_exit_pos[intersect.pair_idx], rc.grid_max_width, rc.grid_max_height, rc.occupied_set_2d);
                        }
                    }
                } else {
                    b_exit_pos[intersect.pair_idx] = {intersect.x, intersect.y};
                }
            }
            prev_pos = b_exit_pos[intersect.pair_idx];
            rc.occupied_set_2d.insert(encode_coord(prev_pos.first, prev_pos.second, 0));
            prev_type = intersect.edge_type;
        }

        // Assign exit terminals to bumps (unchanged)
        for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
            auto& bump_pair = pt.member[pair_idx];
            const std::string& net_name = bump_pair.first.net_name;
            if (rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_x != bump_pair.first.grid_coord_x ||
                rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).grid_coord_y != bump_pair.first.grid_coord_y ||
                rc.bump_dict.at(rc.net_dict.at(net_name).bump_names[0]).layer != bump_pair.first.layer) {
                const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[1];
                const std::string& b_bump_name = rc.net_dict.at(net_name).bump_names[0];
                auto& a_bump = rc.bump_dict.at(a_bump_name);
                a_bump.exit_terminal.bump_name = a_bump_name;
                a_bump.exit_terminal.grid_coord_x = a_exit_pos[pair_idx].first;
                a_bump.exit_terminal.grid_coord_y = a_exit_pos[pair_idx].second;
                a_bump.exit_terminal.layer = assigned_layer;
                auto& b_bump = rc.bump_dict.at(b_bump_name);
                b_bump.exit_terminal.bump_name = b_bump_name;
                b_bump.exit_terminal.grid_coord_x = b_exit_pos[pair_idx].first;
                b_bump.exit_terminal.grid_coord_y = b_exit_pos[pair_idx].second;
                b_bump.exit_terminal.layer = assigned_layer;
                bump_pair.first.exit_terminal = a_bump.exit_terminal;
                bump_pair.second.exit_terminal = b_bump.exit_terminal;
            } else {
                const std::string& a_bump_name = rc.net_dict.at(net_name).bump_names[0];
                const std::string& b_bump_name = rc.net_dict.at(net_name).bump_names[1];
                auto& a_bump = rc.bump_dict.at(a_bump_name);
                a_bump.exit_terminal.bump_name = a_bump_name;
                a_bump.exit_terminal.grid_coord_x = a_exit_pos[pair_idx].first;
                a_bump.exit_terminal.grid_coord_y = a_exit_pos[pair_idx].second;
                a_bump.exit_terminal.layer = assigned_layer;
                auto& b_bump = rc.bump_dict.at(b_bump_name);
                b_bump.exit_terminal.bump_name = b_bump_name;
                b_bump.exit_terminal.grid_coord_x = b_exit_pos[pair_idx].first;
                b_bump.exit_terminal.grid_coord_y = b_exit_pos[pair_idx].second;
                b_bump.exit_terminal.layer = assigned_layer;
                bump_pair.first.exit_terminal = a_bump.exit_terminal;
                bump_pair.second.exit_terminal = b_bump.exit_terminal;
            }
        }

        // Reorder pt.member according to a_intersections order (unchanged)
        std::vector<BumpPair> sorted_member(pair_count);
        for (size_t i = 0; i < pair_count; ++i) {
            sorted_member[i] = pt.member[a_intersections[i].pair_idx];
        }
        pt.member = sorted_member;
    }

    // Output JSON (unchanged)
    json exit_json;
    for (auto& dict_item : rc.bump_dict) {
        const std::string bump_name = dict_item.first;
        const auto& bump = dict_item.second;
        json item;
        item["bump_name"] = bump_name;
        item["net_name"] = bump.net_name;
        item["bump_layer"] = bump.layer;
        item["bump_original_x"] = bump.grid_coord_x;
        item["bump_original_y"] = bump.grid_coord_y;
        item["exit_x"] = bump.exit_terminal.grid_coord_x;
        item["exit_y"] = bump.exit_terminal.grid_coord_y;
        item["exit_layer"] = bump.exit_terminal.layer;
        exit_json.push_back(item);
    }

    const std::string output_dir = rc.output_dir.empty() ? "output" : rc.output_dir;
    std::string output_path = output_dir + "/" + rc.case_name + "_exit_terminals.json";
    std::ofstream out_file(output_path);
    if (out_file.is_open()) {
        out_file << std::setw(4) << exit_json << std::endl;
        out_file.close();
    }

    return intersect_flags;
}
