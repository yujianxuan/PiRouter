#include "ioutils.h"
#include <getopt.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <omp.h>

struct option long_options[] =
{
    {"case", no_argument, 0, 'c'},
    {"thread", no_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {"max_zeros", no_argument, 0, 'z'},
    {"min_rate", no_argument, 0, 'r'},
    {0, 0, 0, 0}
};

const char* ABBRCOMMAND = "htz:r:c:";

const std::string HELP_INFO =
        std::string("Usage: router <case_name> <grid_layout_path> <net_list_path> <output_path> [escape_threads]\n") +
        "Parameters:\n"
        "  <case_name>         Case identifier\n"
        "  <grid_layout_path>  Path to grid_layout input file (required)\n"
        "  <net_list_path>     Path to net_list input file (required)\n"
        "  <output_path>       Path to output file (required)\n"
        "  [escape_threads]    single/serial/false=1, multi/auto/true=team count, or a positive integer\n"
        "Note:\n"
        "  - max_zero is fixed to 2\n"
        "  - default escape thread count is fixed at 13\n";

static std::string get_output_dir_from_path(const std::string& output_path)
{
    size_t pos = output_path.find_last_of("/\\");
    if (pos == std::string::npos) {
        return ".";
    }
    if (pos == 0) {
        return "/";
    }
    return output_path.substr(0, pos);
}

void checkParams(const RunParams& rp)
{
    if (rp.case_name.empty()) {
        IOUtils::terminateProgram("Error: case name is empty!");
    }
    if (rp.grid_layout_path.empty()) {
        IOUtils::terminateProgram("Error: grid_layout path is empty!");
    }
    if (rp.netlist_path.empty()) {
        IOUtils::terminateProgram("Error: netlist path is empty!");
    }
    if (rp.output_path.empty()) {
        IOUtils::terminateProgram("Error: output path is empty!");
    }
}

RunParams parseCmd(int argc, char *argv[])
{
    RunParams rp;
    rp.max_zero = 2;
    rp.escape_threads = 13;

    if (argc < 5 || argc > 6) {
        IOUtils::terminateProgram("Error: incorrect number of parameters!\n" + HELP_INFO);
    }

    rp.case_name = argv[1];
    rp.grid_layout_path = argv[2];
    rp.netlist_path = argv[3];
    rp.output_path = argv[4];
    rp.output_dir = get_output_dir_from_path(rp.output_path);

    if (argc == 6) {
        std::string thread_arg = argv[5];
        if (thread_arg == "true" || thread_arg == "multi" || thread_arg == "auto") {
            rp.escape_threads = 0;
        } else if (thread_arg == "false" || thread_arg == "single" || thread_arg == "serial") {
            rp.escape_threads = 1;
        } else {
            try {
                size_t parsed = 0;
                int requested_threads = std::stoi(thread_arg, &parsed);
                if (parsed != thread_arg.size() || requested_threads <= 0) {
                    IOUtils::terminateProgram("Error: escape thread count must be a positive integer");
                }
                rp.escape_threads = requested_threads;
            } catch (...) {
                IOUtils::terminateProgram("Error: invalid thread argument: " + thread_arg);
            }
        }
    }

    checkParams(rp);
    return rp;
}

void IOUtils::printParams(const RunParams& rp)
{
    IOUtils::print("\n-------------Param Setting-------------");
    IOUtils::print("\nRouting case:\t" + rp.case_name);
    IOUtils::print("\nEscape Threads:\t" + std::to_string(rp.escape_threads));
    IOUtils::print("\nKmeans min rate:\t" + std::to_string(rp.min_rate));
    IOUtils::print("\nKmeans max zero:\t" + std::to_string(rp.max_zero));
    IOUtils::print("\nOutput dir:\t" + rp.output_dir);
    IOUtils::print("\n---------------------------------------\n\n");
}

void IOUtils::terminateProgram(const std::string& errorMessage)
{
    std::cout << errorMessage << std::endl;
    exit(1);
}

void IOUtils::print(const std::string& message)
{
    std::cout << message << std::flush;
}

json coord_to_json(int x, int y, int layer, int layer_count) {
    std::string layer_name = (layer == 0) ? "Bottom" : (layer == layer_count - 1) ? "Top" : "M" + std::to_string(layer);
    return json::array({x, y, layer_name});
}

// ===== 确定性兜底直连（楼梯形折线）=====

// 检查同层一段水平/竖直线 [a..b] 上的格点是否都不被"他 net"占用
// （occupied 里属于 self_cells 的点视为空闲，因同 net 自重叠合法）
static bool seg_free(int x1, int y1, int x2, int y2, int layer,
                     const OccGrid& occ,
                     const std::unordered_set<uint64_t>& self) {
    int dx = (x2 > x1) - (x2 < x1);
    int dy = (y2 > y1) - (y2 < y1);
    int x = x1, y = y1;
    while (true) {
        uint64_t k = encode_coord(x, y, layer);
        if (occ.find(k) != occ.end() && self.find(k) == self.end()) return false;
        if (x == x2 && y == y2) break;
        x += dx; y += dy;
    }
    return true;
}

bool deterministic_connect(
    int sx, int sy, int ex, int ey, int layer,
    int max_width, int max_height,
    OccGrid& occupied_set,
    const std::unordered_set<uint64_t>& self_cells,
    std::vector<Seg>& segs_out
) {
    (void)max_width;
    // 楼梯形：sx,sy --竖直--> sx,gy --水平--> ex,gy --竖直--> ex,ey
    // 找一条 gy，使三段都不撞他 net。从 sy/ey 邻近向外扩。
    int lo = std::min(sy, ey), hi = std::max(sy, ey);
    for (int d = 0; d < max_height; ++d) {
        for (int gy : { hi + d, lo - d }) {
            if (gy < 0 || gy >= max_height) continue;
            if (!seg_free(sx, sy, sx, gy, layer, occupied_set, self_cells)) continue;
            if (!seg_free(sx, gy, ex, gy, layer, occupied_set, self_cells)) continue;
            if (!seg_free(ex, gy, ex, ey, layer, occupied_set, self_cells)) continue;
            auto emit = [&](int x1, int y1, int x2, int y2) {
                if (x1 == x2 && y1 == y2) return;
                segs_out.push_back(Seg{ ResultPoint{x1, y1, layer}, ResultPoint{x2, y2, layer} });
                int ddx = (x2 > x1) - (x2 < x1), ddy = (y2 > y1) - (y2 < y1);
                int x = x1, y = y1;
                while (true) {
                    occupied_set.insert(encode_coord(x, y, layer));
                    if (x == x2 && y == y2) break;
                    x += ddx; y += ddy;
                }
            };
            emit(sx, sy, sx, gy);
            emit(sx, gy, ex, gy);
            emit(ex, gy, ex, ey);
            return true;
        }
    }
    return false;
}

// 换层 U 形兜底：两 bump 同 xy 之一时，via 到中间空层 L'，在 L' 上走 U 形绕到空闲列/行再回来。
// 比同层楼梯更强：在所有中间层 × 多个偏移上找完全空闲通道。成功写 segs_out + occupied_set。
// bump 在 bump_layer 层，中间层 L' 取 1..layer_count-2（避开 Bottom=0 拥挤、Top=layer_count-1）。
bool layered_u_connect(
    int sx, int sy, int ex, int ey, int bump_layer, int layer_count,
    int max_width, int max_height,
    OccGrid& occupied_set,
    const std::unordered_set<uint64_t>& self_cells,
    std::vector<Seg>& segs_out
) {
    auto via_free = [&](int x, int y, int l1, int l2) {
        int lo = std::min(l1,l2), hi = std::max(l1,l2);
        for (int l = lo; l <= hi; ++l) {
            uint64_t k = encode_coord(x, y, l);
            if (occupied_set.find(k) != occupied_set.end() && self_cells.find(k) == self_cells.end()) return false;
        }
        return true;
    };
    bool same_x = (sx == ex);
    bool same_y = (sy == ey);
    if (!same_x && !same_y) return false;  // 只处理同 x 或同 y 的规整结构

    for (int L = 1; L < layer_count - 1; ++L) {       // 中间层
        if (L == bump_layer) continue;
        // via 列必须能从 bump_layer 通到 L
        if (!via_free(sx, sy, bump_layer, L) || !via_free(ex, ey, bump_layer, L)) continue;
        for (int off = 5; off <= 400; off += 5) {     // 偏移找空闲通道
            for (int sign : {-1, 1}) {
                if (same_x) {
                    int gx = sx + sign * off;
                    if (gx < 0 || gx >= max_width) continue;
                    // U形: (sx,sy)->(gx,sy)->(gx,ey)->(ex,ey) 全在 L 层
                    if (!seg_free(sx, sy, gx, sy, L, occupied_set, self_cells)) continue;
                    if (!seg_free(gx, sy, gx, ey, L, occupied_set, self_cells)) continue;
                    if (!seg_free(gx, ey, ex, ey, L, occupied_set, self_cells)) continue;
                    auto emitL = [&](int x1,int y1,int x2,int y2){
                        if (x1==x2&&y1==y2) return;
                        segs_out.push_back(Seg{ResultPoint{x1,y1,L},ResultPoint{x2,y2,L}});
                        int ddx=(x2>x1)-(x2<x1),ddy=(y2>y1)-(y2<y1),x=x1,y=y1;
                        while(true){occupied_set.insert(encode_coord(x,y,L));if(x==x2&&y==y2)break;x+=ddx;y+=ddy;}
                    };
                    auto emitVia = [&](int x,int y){
                        int lo=std::min(bump_layer,L),hi=std::max(bump_layer,L);
                        segs_out.push_back(Seg{ResultPoint{x,y,bump_layer},ResultPoint{x,y,L}});
                        for(int l=lo;l<=hi;++l) occupied_set.insert(encode_coord(x,y,l));
                    };
                    emitVia(sx, sy);
                    emitL(sx, sy, gx, sy);
                    emitL(gx, sy, gx, ey);
                    emitL(gx, ey, ex, ey);
                    emitVia(ex, ey);
                    return true;
                } else { // same_y
                    int gy = sy + sign * off;
                    if (gy < 0 || gy >= max_height) continue;
                    if (!seg_free(sx, sy, sx, gy, L, occupied_set, self_cells)) continue;
                    if (!seg_free(sx, gy, ex, gy, L, occupied_set, self_cells)) continue;
                    if (!seg_free(ex, gy, ex, ey, L, occupied_set, self_cells)) continue;
                    auto emitL = [&](int x1,int y1,int x2,int y2){
                        if (x1==x2&&y1==y2) return;
                        segs_out.push_back(Seg{ResultPoint{x1,y1,L},ResultPoint{x2,y2,L}});
                        int ddx=(x2>x1)-(x2<x1),ddy=(y2>y1)-(y2<y1),x=x1,y=y1;
                        while(true){occupied_set.insert(encode_coord(x,y,L));if(x==x2&&y==y2)break;x+=ddx;y+=ddy;}
                    };
                    auto emitVia = [&](int x,int y){
                        int lo=std::min(bump_layer,L),hi=std::max(bump_layer,L);
                        segs_out.push_back(Seg{ResultPoint{x,y,bump_layer},ResultPoint{x,y,L}});
                        for(int l=lo;l<=hi;++l) occupied_set.insert(encode_coord(x,y,l));
                    };
                    emitVia(sx, sy);
                    emitL(sx, sy, sx, gy);
                    emitL(sx, gy, ex, gy);
                    emitL(ex, gy, ex, ey);
                    emitVia(ex, ey);
                    return true;
                }
            }
        }
    }
    return false;
}

bool layered_stair_connect(
    int sx, int sy, int ex, int ey, int bump_layer, int layer_count,
    int max_width, int max_height,
    OccGrid& occupied_set,
    const std::unordered_set<uint64_t>& self_cells,
    std::vector<Seg>& segs_out
) {
    (void)max_width;
    auto via_free = [&](int x, int y, int l1, int l2) {
        uint64_t k1 = encode_coord(x, y, l1);
        uint64_t k2 = encode_coord(x, y, l2);
        if (occupied_set.find(k1) != occupied_set.end() && self_cells.find(k1) == self_cells.end()) {
            return false;
        }
        if (occupied_set.find(k2) != occupied_set.end() && self_cells.find(k2) == self_cells.end()) {
            return false;
        }
        return true;
    };

    auto emit_line = [&](int x1, int y1, int x2, int y2, int layer) {
        if (x1 == x2 && y1 == y2) return;
        segs_out.push_back(Seg{ResultPoint{x1, y1, layer}, ResultPoint{x2, y2, layer}});
        int dx = (x2 > x1) - (x2 < x1);
        int dy = (y2 > y1) - (y2 < y1);
        int x = x1;
        int y = y1;
        while (true) {
            occupied_set.insert(encode_coord(x, y, layer));
            if (x == x2 && y == y2) break;
            x += dx;
            y += dy;
        }
    };

    auto emit_via = [&](int x, int y, int from_layer, int to_layer) {
        if (from_layer == to_layer) return;
        segs_out.push_back(Seg{ResultPoint{x, y, from_layer}, ResultPoint{x, y, to_layer}});
        occupied_set.insert(encode_coord(x, y, from_layer));
        occupied_set.insert(encode_coord(x, y, to_layer));
    };

    for (int L = 1; L < layer_count - 1; ++L) {
        if (L == bump_layer) continue;
        if (!via_free(sx, sy, bump_layer, L) || !via_free(ex, ey, bump_layer, L)) continue;
        int lo = std::min(sy, ey);
        int hi = std::max(sy, ey);
        for (int d = 0; d < max_height; ++d) {
            for (int gy : {hi + d, lo - d}) {
                if (gy < 0 || gy >= max_height) continue;
                if (!seg_free(sx, sy, sx, gy, L, occupied_set, self_cells)) continue;
                if (!seg_free(sx, gy, ex, gy, L, occupied_set, self_cells)) continue;
                if (!seg_free(ex, gy, ex, ey, L, occupied_set, self_cells)) continue;
                emit_via(sx, sy, bump_layer, L);
                emit_line(sx, sy, sx, gy, L);
                emit_line(sx, gy, ex, gy, L);
                emit_line(ex, gy, ex, ey, L);
                emit_via(ex, ey, L, bump_layer);
                return true;
            }
        }
    }
    return false;
}

// 移植 evaluator.py 的内角公式：返回 0~180 度，0=反向，180=同向
static double angle_between_pts(const ResultPoint& a, const ResultPoint& b, const ResultPoint& c) {
    double v1x = a.x - b.x, v1y = a.y - b.y;
    double v2x = c.x - b.x, v2y = c.y - b.y;
    if ((std::abs(v1x) < 1e-6 && std::abs(v1y) < 1e-6) ||
        (std::abs(v2x) < 1e-6 && std::abs(v2y) < 1e-6))
        return 180.0;
    double dot = v1x * v2x + v1y * v2y;
    double crs = v1x * v2y - v1y * v2x;
    return std::atan2(std::abs(crs), dot) * 180.0 / M_PI;
}

static bool net_contains_point(
    const std::vector<Seg>& segs, int x, int y, int layer
) {
    for (const Seg& seg : segs) {
        const ResultPoint& a = seg.first;
        const ResultPoint& b = seg.second;
        if (a.layer != b.layer) {
            if ((x == a.x && y == a.y && layer == a.layer) ||
                (x == b.x && y == b.y && layer == b.layer)) return true;
            continue;
        }
        if (layer != a.layer || x < std::min(a.x, b.x) || x > std::max(a.x, b.x) ||
            y < std::min(a.y, b.y) || y > std::max(a.y, b.y)) continue;
        int dx = b.x - a.x;
        int dy = b.y - a.y;
        if ((dx == 0 && x == a.x) || (dy == 0 && y == a.y) ||
            (std::abs(dx) == std::abs(dy) && std::abs(x - a.x) == std::abs(y - a.y))) {
            return true;
        }
    }
    return false;
}

// A destination may be free or already owned by this net.
static bool via_dest_free(
    int x, int y, int layer,
    const OccGrid& occupied_set,
    const std::unordered_set<uint64_t>& self_points,
    const std::vector<Seg>& self_segments
) {
    uint64_t key = encode_coord(x, y, layer);
    return occupied_set.find(key) == occupied_set.end() ||
           self_points.find(key) != self_points.end() ||
           net_contains_point(self_segments, x, y, layer);
}

static double segment_len(const ResultPoint& a, const ResultPoint& b) {
    if (a.x == b.x && a.y == b.y && a.layer != b.layer) return 0.0;
    return std::hypot(static_cast<double>(b.x - a.x),
                      static_cast<double>(b.y - a.y));
}

static bool is_queen_line(const ResultPoint& a, const ResultPoint& b) {
    if (a.layer != b.layer) return a.x == b.x && a.y == b.y;
    int dx = std::abs(b.x - a.x);
    int dy = std::abs(b.y - a.y);
    return dx == 0 || dy == 0 || dx == dy;
}

static void add_segment_occupancy_keys(
    const ResultPoint& a,
    const ResultPoint& b,
    std::unordered_set<uint64_t>& keys
) {
    if (a.x == b.x && a.y == b.y && a.layer != b.layer) {
        keys.insert(encode_coord(a.x, a.y, a.layer));
        keys.insert(encode_coord(b.x, b.y, b.layer));
        return;
    }
    if (a.layer != b.layer) return;
    int dx = b.x - a.x;
    int dy = b.y - a.y;
    int steps = std::max(std::abs(dx), std::abs(dy));
    if (steps == 0) {
        keys.insert(encode_coord(a.x, a.y, a.layer));
        return;
    }
    int step_x = (dx > 0) - (dx < 0);
    int step_y = (dy > 0) - (dy < 0);
    int x = a.x;
    int y = a.y;
    for (int i = 0; i <= steps; ++i) {
        keys.insert(encode_coord(x, y, a.layer));
        if (i < steps && step_x != 0 && step_y != 0) {
            keys.insert(encode_coord(x + step_x * 0.5f,
                                     y + step_y * 0.5f,
                                     a.layer));
        }
        x += step_x;
        y += step_y;
    }
}

static void replace_path_occupancy(
    const std::unordered_set<uint64_t>& stale_keys,
    const std::unordered_set<uint64_t>& final_keys,
    OccGrid& occupied_set
) {
    for (uint64_t key : stale_keys) {
        occupied_set.erase(key);
    }
    for (uint64_t key : final_keys) {
        occupied_set.insert(key);
    }
}

static void replace_path_occupancy(
    const std::unordered_set<uint64_t>& stale_keys,
    const std::vector<uint64_t>& final_keys,
    OccGrid& occupied_set
) {
    for (uint64_t key : stale_keys) {
        occupied_set.erase(key);
    }
    for (uint64_t key : final_keys) {
        occupied_set.insert(key);
    }
}

static void insert_segment_occupancy(
    const ResultPoint& a,
    const ResultPoint& b,
    OccGrid& occupied_set
) {
    if (a.x == b.x && a.y == b.y && a.layer != b.layer) {
        occupied_set.insert(encode_coord(a.x, a.y, a.layer));
        occupied_set.insert(encode_coord(b.x, b.y, b.layer));
        return;
    }
    if (a.layer != b.layer) return;
    const int dx = b.x - a.x;
    const int dy = b.y - a.y;
    const int steps = std::max(std::abs(dx), std::abs(dy));
    const int step_x = (dx > 0) - (dx < 0);
    const int step_y = (dy > 0) - (dy < 0);
    int x = a.x;
    int y = a.y;
    for (int step = 0; step <= steps; ++step) {
        occupied_set.insert(encode_coord(x, y, a.layer));
        if (step < steps && step_x != 0 && step_y != 0) {
            occupied_set.insert(encode_coord(
                x + step_x * 0.5f, y + step_y * 0.5f, a.layer));
        }
        x += step_x;
        y += step_y;
    }
}

static void collect_path_occupancy_keys(
    const std::vector<Seg>& segs,
    std::unordered_set<uint64_t>& keys
) {
    keys.clear();
    for (const Seg& seg : segs) {
        add_segment_occupancy_keys(seg.first, seg.second, keys);
    }
}

static bool shortcut_clear(
    const ResultPoint& a,
    const ResultPoint& b,
    const OccGrid& occupied_set,
    const std::unordered_set<uint64_t>& self_keys
) {
    if (a.layer != b.layer || !is_queen_line(a, b)) return false;
    int dx = b.x - a.x;
    int dy = b.y - a.y;
    int steps = std::max(std::abs(dx), std::abs(dy));
    if (steps == 0) return false;
    int step_x = (dx > 0) - (dx < 0);
    int step_y = (dy > 0) - (dy < 0);
    int x = a.x;
    int y = a.y;
    for (int i = 0; i <= steps; ++i) {
        uint64_t key = encode_coord(x, y, a.layer);
        if (occupied_set.find(key) != occupied_set.end() &&
            self_keys.find(key) == self_keys.end()) {
            return false;
        }
        x += step_x;
        y += step_y;
    }
    return true;
}

static bool boundary_angle_ok(
    const std::vector<ResultPoint>& vertices,
    size_t from,
    size_t to
) {
    const ResultPoint& a = vertices[from];
    const ResultPoint& b = vertices[to];
    if (from > 0) {
        const ResultPoint& prev = vertices[from - 1];
        if (prev.layer == a.layer && a.layer == b.layer &&
            (prev.x != a.x || prev.y != a.y)) {
            double ang = angle_between_pts(prev, a, b);
            if (ang > 0.0 && ang < 90.0) return false;
        }
    }
    if (to + 1 < vertices.size()) {
        const ResultPoint& next = vertices[to + 1];
        if (a.layer == b.layer && b.layer == next.layer &&
            (next.x != b.x || next.y != b.y)) {
            double ang = angle_between_pts(a, b, next);
            if (ang > 0.0 && ang < 90.0) return false;
        }
    }
    return true;
}

static bool same_direction(const ResultPoint& a, const ResultPoint& b,
                           const ResultPoint& c) {
    if (a.layer != b.layer || b.layer != c.layer) return false;
    int dx1 = (b.x > a.x) - (b.x < a.x);
    int dy1 = (b.y > a.y) - (b.y < a.y);
    int dx2 = (c.x > b.x) - (c.x < b.x);
    int dy2 = (c.y > b.y) - (c.y < b.y);
    return dx1 == dx2 && dy1 == dy2;
}

std::unordered_set<std::string> shorten_routing_paths(
    std::unordered_map<std::string, std::vector<Seg>>& routing_paths,
    OccGrid& occupied_set,
    const std::unordered_set<std::string>& skip_nets
) {
    size_t changed_nets = 0;
    std::unordered_set<std::string> changed_net_names;
    size_t shortcut_count = 0;
    double saved_length = 0.0;
    constexpr size_t MAX_JUMP_VERTICES = 64;
    constexpr double MIN_GAIN = 0.25;

    struct Edge {
        int a;
        int b;
        double len;
    };

    for (auto& kv : routing_paths) {
        if (skip_nets.count(kv.first)) continue;
        std::vector<Seg>& segs = kv.second;
        if (segs.size() < 3) continue;

        std::unordered_set<uint64_t> self_keys;
        self_keys.reserve(segs.size() * 8);
        for (const Seg& seg : segs) {
            add_segment_occupancy_keys(seg.first, seg.second, self_keys);
        }

        std::vector<ResultPoint> points;
        std::unordered_map<uint64_t, int> point_index;
        std::vector<Edge> edges;
        std::vector<std::vector<std::pair<int, int>>> adj;

        auto get_index = [&](const ResultPoint& p) {
            uint64_t key = encode_coord(p.x, p.y, p.layer);
            auto it = point_index.find(key);
            if (it != point_index.end()) return it->second;
            int idx = static_cast<int>(points.size());
            point_index.emplace(key, idx);
            points.push_back(p);
            adj.emplace_back();
            return idx;
        };

        bool valid_graph = true;
        for (const Seg& seg : segs) {
            if (!is_queen_line(seg.first, seg.second)) {
                valid_graph = false;
                break;
            }
            int a = get_index(seg.first);
            int b = get_index(seg.second);
            if (a == b) continue;
            int ei = static_cast<int>(edges.size());
            edges.push_back(Edge{a, b, segment_len(seg.first, seg.second)});
            adj[a].push_back({b, ei});
            adj[b].push_back({a, ei});
            if (adj[a].size() > 2 || adj[b].size() > 2) {
                valid_graph = false;
                break;
            }
        }
        if (!valid_graph || edges.size() < 3) continue;

        std::vector<int> endpoints;
        endpoints.reserve(2);
        for (size_t i = 0; i < adj.size(); ++i) {
            if (adj[i].empty()) continue;
            if (adj[i].size() == 1) endpoints.push_back(static_cast<int>(i));
            if (adj[i].size() > 2) {
                valid_graph = false;
                break;
            }
        }
        if (!valid_graph || endpoints.size() != 2) continue;

        std::vector<ResultPoint> vertices;
        std::vector<double> prefix_len;
        vertices.reserve(edges.size() + 1);
        prefix_len.reserve(edges.size() + 1);
        std::vector<char> used_edge(edges.size(), 0);
        int prev = -1;
        int cur = endpoints[0];
        vertices.push_back(points[cur]);
        prefix_len.push_back(0.0);
        while (cur != endpoints[1]) {
            int next = -1;
            int edge_idx = -1;
            for (const auto& item : adj[cur]) {
                if (item.first == prev) continue;
                next = item.first;
                edge_idx = item.second;
                break;
            }
            if (next < 0 || edge_idx < 0 || used_edge[edge_idx]) {
                valid_graph = false;
                break;
            }
            used_edge[edge_idx] = 1;
            prefix_len.push_back(prefix_len.back() + edges[edge_idx].len);
            vertices.push_back(points[next]);
            prev = cur;
            cur = next;
        }
        if (!valid_graph || vertices.size() != edges.size() + 1) continue;

        bool changed = false;
        std::vector<ResultPoint> new_vertices;
        new_vertices.reserve(vertices.size());
        size_t i = 0;
        new_vertices.push_back(vertices[0]);
        while (i + 1 < vertices.size()) {
            size_t best_j = i + 1;
            double best_gain = 0.0;
            size_t limit = std::min(vertices.size() - 1, i + MAX_JUMP_VERTICES);
            for (size_t j = i + 2; j <= limit; ++j) {
                const ResultPoint& a = vertices[i];
                const ResultPoint& b = vertices[j];
                if (a.layer != b.layer) continue;
                if (!is_queen_line(a, b)) continue;
                double old_len = prefix_len[j] - prefix_len[i];
                double new_len = segment_len(a, b);
                double gain = old_len - new_len;
                if (gain <= best_gain + MIN_GAIN) continue;
                if (!boundary_angle_ok(vertices, i, j)) continue;
                if (!shortcut_clear(a, b, occupied_set, self_keys)) continue;
                best_gain = gain;
                best_j = j;
            }
            if (best_j != i + 1) {
                changed = true;
                shortcut_count++;
                saved_length += best_gain;
            }
            new_vertices.push_back(vertices[best_j]);
            i = best_j;
        }
        if (!changed) continue;

        std::vector<ResultPoint> compressed;
        compressed.reserve(new_vertices.size());
        for (const ResultPoint& p : new_vertices) {
            if (compressed.empty()) {
                compressed.push_back(p);
                continue;
            }
            if (p.x == compressed.back().x &&
                p.y == compressed.back().y &&
                p.layer == compressed.back().layer) {
                continue;
            }
            if (compressed.size() >= 2 &&
                same_direction(compressed[compressed.size() - 2],
                               compressed.back(), p)) {
                compressed.back() = p;
            } else {
                compressed.push_back(p);
            }
        }
        if (compressed.size() < 2) continue;

        std::vector<Seg> replacement;
        std::unordered_set<uint64_t> replacement_keys;
        replacement_keys.reserve(self_keys.size());
        replacement.reserve(compressed.size() - 1);
        for (size_t k = 1; k < compressed.size(); ++k) {
            if (compressed[k - 1].x == compressed[k].x &&
                compressed[k - 1].y == compressed[k].y &&
                compressed[k - 1].layer == compressed[k].layer) {
                continue;
            }
            replacement.emplace_back(compressed[k - 1], compressed[k]);
            add_segment_occupancy_keys(compressed[k - 1], compressed[k],
                                       replacement_keys);
        }
        if (replacement.empty()) continue;
        segs.swap(replacement);
        replace_path_occupancy(self_keys, replacement_keys, occupied_set);
        changed_nets++;
        changed_net_names.insert(kv.first);
    }

    std::cout << "[Shorten] changed " << changed_nets
              << " nets, shortcuts " << shortcut_count
              << ", estimated saved length " << std::fixed << std::setprecision(3)
              << saved_length << std::endl;
    return changed_net_names;
}

static int segment_via_count(const Seg& seg) {
    const ResultPoint& a = seg.first;
    const ResultPoint& b = seg.second;
    if (a.layer == b.layer) return 0;
    if (a.x == b.x && a.y == b.y) return std::abs(a.layer - b.layer);
    return 0;
}

static long long net_via_count(const std::vector<Seg>& segs) {
    long long total = 0;
    for (const Seg& seg : segs) {
        total += segment_via_count(seg);
    }
    return total;
}

static void rebuild_occupied_from_paths(
    const RoutingCase& rc,
    const std::unordered_map<std::string, std::vector<Seg>>& routing_paths,
    OccGrid& occupied_set
) {
    occupied_set.clear_blocks();
    int max_coord = std::max(rc.grid_max_width, rc.grid_max_height) + 10000;
    occupied_set.init(max_coord, rc.layer_count);
    for (const std::string& net_name : rc.net_name_list) {
        auto it = routing_paths.find(net_name);
        if (it == routing_paths.end()) continue;
        for (const Seg& seg : it->second) {
            insert_segment_occupancy(
                seg.first, seg.second, occupied_set);
        }
    }
}

static bool candidate_segments_clear_vector(
    const std::vector<Seg>& candidate,
    const OccGrid& occupied_set,
    const std::unordered_set<uint64_t>& self_keys,
    std::vector<uint64_t>& candidate_keys
) {
    candidate_keys.clear();
    auto append_key = [&](uint64_t key) {
        if (occupied_set.find(key) != occupied_set.end() &&
            self_keys.find(key) == self_keys.end()) {
            return false;
        }
        candidate_keys.push_back(key);
        return true;
    };
    for (const Seg& seg : candidate) {
        const ResultPoint& a = seg.first;
        const ResultPoint& b = seg.second;
        if (!is_queen_line(a, b)) return false;
        if (a.x == b.x && a.y == b.y && a.layer != b.layer) {
            if (!append_key(encode_coord(a.x, a.y, a.layer)) ||
                !append_key(encode_coord(b.x, b.y, b.layer))) {
                return false;
            }
            continue;
        }
        if (a.layer != b.layer) return false;
        const int dx = b.x - a.x;
        const int dy = b.y - a.y;
        const int steps = std::max(std::abs(dx), std::abs(dy));
        const int step_x = (dx > 0) - (dx < 0);
        const int step_y = (dy > 0) - (dy < 0);
        int x = a.x;
        int y = a.y;
        for (int step = 0; step <= steps; ++step) {
            if (!append_key(encode_coord(x, y, a.layer))) return false;
            if (step < steps && step_x != 0 && step_y != 0 &&
                !append_key(encode_coord(
                    x + step_x * 0.5f,
                    y + step_y * 0.5f,
                    a.layer))) {
                return false;
            }
            x += step_x;
            y += step_y;
        }
    }
    return true;
}

static bool build_layer_shift_candidate(
    const RoutingCase& rc,
    const std::string& net_name,
    const std::vector<Seg>& segs,
    int target_layer,
    std::vector<Seg>& candidate
) {
    candidate.clear();
    candidate.reserve(segs.size() + 2);
    for (const Seg& seg : segs) {
        ResultPoint a{seg.first.x, seg.first.y, target_layer};
        ResultPoint b{seg.second.x, seg.second.y, target_layer};
        if (a.x == b.x && a.y == b.y) continue;
        if (!is_queen_line(a, b)) return false;
        candidate.emplace_back(a, b);
    }

    const Net& net = rc.net_dict.at(net_name);
    for (const std::string& bump_name : net.bump_names) {
        const Bump& bump = rc.bump_dict.at(bump_name);
        if (bump.layer != target_layer) {
            candidate.emplace_back(
                ResultPoint{bump.grid_coord_x, bump.grid_coord_y, bump.layer},
                ResultPoint{bump.grid_coord_x, bump.grid_coord_y, target_layer});
        }
    }
    return !candidate.empty();
}

void reduce_pin_stub_vias(
    const RoutingCase& rc,
    std::unordered_map<std::string, std::vector<Seg>>& routing_paths,
    OccGrid& occupied_set
) {
    const int top_layer = rc.layer_count - 1;
    constexpr long long MIN_ELIGIBLE_VIAS_FOR_LAYER_SHIFT = 5000;

    struct ViaCandidate {
        std::string net_name;
        long long vias;
        int pin_layer;
    };

    std::vector<ViaCandidate> ordered_nets;
    long long eligible_vias = 0;
    for (const std::string& net_name : rc.net_name_list) {
        auto path_it = routing_paths.find(net_name);
        if (path_it == routing_paths.end() || path_it->second.empty()) continue;
        auto net_it = rc.net_dict.find(net_name);
        if (net_it == rc.net_dict.end()) continue;

        const Net& net = net_it->second;
        if (net.bump_names.empty()) continue;
        int pin_layer = rc.bump_dict.at(net.bump_names.front()).layer;
        bool pure_pin_layer = (pin_layer == 0 || pin_layer == top_layer);
        for (const std::string& bump_name : net.bump_names) {
            auto bump_it = rc.bump_dict.find(bump_name);
            if (bump_it == rc.bump_dict.end() ||
                bump_it->second.layer != pin_layer) {
                pure_pin_layer = false;
                break;
            }
        }
        if (!pure_pin_layer) continue;

        long long vias = net_via_count(path_it->second);
        if (vias == 0) continue;
        eligible_vias += vias;
        ordered_nets.push_back(ViaCandidate{net_name, vias, pin_layer});
    }

    rebuild_occupied_from_paths(rc, routing_paths, occupied_set);

    if (eligible_vias < MIN_ELIGIBLE_VIAS_FOR_LAYER_SHIFT) {
        std::cout << "[ViaReduce] skipped, eligible vias " << eligible_vias
                  << " below threshold "
                  << MIN_ELIGIBLE_VIAS_FOR_LAYER_SHIFT << std::endl;
        return;
    }

    size_t changed_nets = 0;
    long long saved_vias = 0;
    std::sort(ordered_nets.begin(), ordered_nets.end(),
              [](const ViaCandidate& a, const ViaCandidate& b) {
        if (a.vias != b.vias) return a.vias > b.vias;
        return a.net_name < b.net_name;
    });

    std::unordered_set<uint64_t> self_keys;
    std::vector<uint64_t> candidate_keys;
    std::vector<Seg> candidate;

    for (const ViaCandidate& via_candidate : ordered_nets) {
        const std::string& net_name = via_candidate.net_name;
        auto path_it = routing_paths.find(net_name);
        if (path_it == routing_paths.end()) continue;
        std::vector<Seg>& segs = path_it->second;
        if (segs.empty()) continue;

        const Net& net = rc.net_dict.at(net_name);
        int pin_layer = via_candidate.pin_layer;
        bool pure_pin_layer = (pin_layer == 0 || pin_layer == top_layer);
        for (const std::string& bump_name : net.bump_names) {
            auto bump_it = rc.bump_dict.find(bump_name);
            if (bump_it == rc.bump_dict.end() ||
                bump_it->second.layer != pin_layer) {
                pure_pin_layer = false;
                break;
            }
        }
        if (!pure_pin_layer) continue;

        long long old_vias = net_via_count(segs);
        if (old_vias == 0) continue;

        self_keys.clear();
        for (const Seg& seg : segs) {
            add_segment_occupancy_keys(seg.first, seg.second, self_keys);
        }

        std::vector<int> target_layers;
        target_layers.reserve(rc.layer_count);
        if (pin_layer == top_layer) {
            for (int layer = top_layer; layer >= 0; --layer) {
                target_layers.push_back(layer);
            }
        } else {
            for (int layer = 0; layer < rc.layer_count; ++layer) {
                target_layers.push_back(layer);
            }
        }

        for (int target_layer : target_layers) {
            if (!build_layer_shift_candidate(rc, net_name, segs, target_layer, candidate)) continue;
            long long new_vias = net_via_count(candidate);
            if (new_vias >= old_vias) continue;
            if (!candidate_segments_clear_vector(
                    candidate, occupied_set, self_keys, candidate_keys)) {
                continue;
            }

            segs.swap(candidate);
            replace_path_occupancy(self_keys, candidate_keys, occupied_set);
            changed_nets++;
            saved_vias += old_vias - new_vias;
            break;
        }
    }

    std::cout << "[ViaReduce] changed " << changed_nets
              << " nets, estimated saved vias " << saved_vias << std::endl;
}

static bool planar_segment(const Seg& seg) {
    return seg.first.layer == seg.second.layer &&
           (seg.first.x != seg.second.x || seg.first.y != seg.second.y);
}

static bool via_between_at(
    const Seg& seg,
    const ResultPoint& p,
    int layer,
    int& other_layer
) {
    const ResultPoint& a = seg.first;
    const ResultPoint& b = seg.second;
    if (a.x != b.x || a.y != b.y || a.layer == b.layer) return false;
    if (a.x != p.x || a.y != p.y) return false;
    if (a.layer == layer) {
        other_layer = b.layer;
        return true;
    }
    if (b.layer == layer) {
        other_layer = a.layer;
        return true;
    }
    return false;
}

static long long two_boundary_vias(
    int run_layer,
    int before_layer,
    int after_layer
) {
    return std::abs(run_layer - before_layer) +
           std::abs(run_layer - after_layer);
}

static bool same_point(const ResultPoint& a, const ResultPoint& b) {
    return a.x == b.x && a.y == b.y && a.layer == b.layer;
}

long long reduce_internal_run_vias(
    const RoutingCase& rc,
    std::unordered_map<std::string, std::vector<Seg>>& routing_paths,
    OccGrid& occupied_set
) {
    constexpr long long MIN_TOTAL_VIAS_FOR_RUN_SHIFT = 50000;

    long long total_vias = 0;
    for (const auto& kv : routing_paths) {
        total_vias += net_via_count(kv.second);
    }
    if (total_vias < MIN_TOTAL_VIAS_FOR_RUN_SHIFT) {
        std::cout << "[RunViaReduce] skipped, total vias " << total_vias
                  << " below threshold "
                  << MIN_TOTAL_VIAS_FOR_RUN_SHIFT << std::endl;
        return 0;
    }

    size_t changed_runs = 0;
    size_t changed_nets = 0;
    long long saved_vias = 0;

    std::unordered_set<uint64_t> self_keys;
    std::vector<uint64_t> candidate_keys;
    std::unordered_set<uint64_t> final_keys;
    std::vector<Seg> candidate_run;
    std::vector<Seg> candidate_replacement;

    for (auto& kv : routing_paths) {
        std::vector<Seg>& segs = kv.second;
        if (segs.size() < 3) continue;

        bool net_changed = false;
        collect_path_occupancy_keys(segs, self_keys);

        size_t i = 1;
        while (i + 1 < segs.size()) {
            if (!planar_segment(segs[i])) {
                ++i;
                continue;
            }

            const int run_layer = segs[i].first.layer;
            const ResultPoint run_start = segs[i].first;
            ResultPoint run_end = segs[i].second;
            size_t j = i + 1;
            bool contiguous = true;
            while (j + 1 < segs.size() && planar_segment(segs[j])) {
                if (segs[j].first.layer != run_layer ||
                    !same_point(segs[j].first, run_end)) {
                    contiguous = false;
                    break;
                }
                run_end = segs[j].second;
                ++j;
            }
            if (!contiguous) {
                ++i;
                continue;
            }

            int before_layer = -1;
            int after_layer = -1;
            if (!via_between_at(segs[i - 1], run_start, run_layer, before_layer) ||
                !via_between_at(segs[j], run_end, run_layer, after_layer)) {
                i = j;
                continue;
            }

            long long old_vias = two_boundary_vias(run_layer, before_layer, after_layer);
            std::vector<std::pair<long long, int>> layer_candidates;
            layer_candidates.reserve(rc.layer_count - 1);
            for (int layer = 0; layer < rc.layer_count; ++layer) {
                if (layer == run_layer) continue;
                long long new_vias = two_boundary_vias(layer, before_layer, after_layer);
                if (new_vias < old_vias) {
                    layer_candidates.emplace_back(old_vias - new_vias, layer);
                }
            }
            if (layer_candidates.empty()) {
                i = j;
                continue;
            }
            std::sort(layer_candidates.begin(), layer_candidates.end(),
                      [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
            });

            for (const auto& layer_candidate : layer_candidates) {
                long long via_gain = layer_candidate.first;
                int target_layer = layer_candidate.second;

                candidate_run.clear();
                candidate_run.reserve(j - i);
                for (size_t k = i; k < j; ++k) {
                    candidate_run.emplace_back(
                        ResultPoint{segs[k].first.x, segs[k].first.y, target_layer},
                        ResultPoint{segs[k].second.x, segs[k].second.y, target_layer});
                }
                Seg before_replacement;
                if (same_point(segs[i - 1].second, run_start)) {
                    before_replacement = Seg{
                        ResultPoint{run_start.x, run_start.y, before_layer},
                        ResultPoint{run_start.x, run_start.y, target_layer}};
                } else {
                    before_replacement = Seg{
                        ResultPoint{run_start.x, run_start.y, target_layer},
                        ResultPoint{run_start.x, run_start.y, before_layer}};
                }
                Seg after_replacement;
                if (same_point(segs[j].first, run_end)) {
                    after_replacement = Seg{
                        ResultPoint{run_end.x, run_end.y, target_layer},
                        ResultPoint{run_end.x, run_end.y, after_layer}};
                } else {
                    after_replacement = Seg{
                        ResultPoint{run_end.x, run_end.y, after_layer},
                        ResultPoint{run_end.x, run_end.y, target_layer}};
                }

                // This transformation preserves graph connectivity locally:
                // the two rewritten vias retain their outer endpoints, and
                // connect them to the unchanged XY run on target_layer. Check
                // the complete replacement (including both vias) for legality
                // and occupancy, then update only the affected segment range.
                candidate_replacement.clear();
                candidate_replacement.reserve(candidate_run.size() + 2);
                candidate_replacement.push_back(before_replacement);
                candidate_replacement.insert(
                    candidate_replacement.end(),
                    candidate_run.begin(), candidate_run.end());
                candidate_replacement.push_back(after_replacement);
                if (!candidate_segments_clear_vector(
                        candidate_replacement, occupied_set,
                        self_keys, candidate_keys)) {
                    continue;
                }

                segs[i - 1] = before_replacement;
                for (size_t k = i; k < j; ++k) {
                    segs[k] = candidate_run[k - i];
                }
                segs[j] = after_replacement;
                for (uint64_t key : candidate_keys) {
                    occupied_set.insert(key);
                    self_keys.insert(key);
                }
                self_keys.insert(encode_coord(run_start.x, run_start.y, target_layer));
                self_keys.insert(encode_coord(run_end.x, run_end.y, target_layer));

                changed_runs++;
                saved_vias += via_gain;
                net_changed = true;
                break;
            }
            i = j;
        }
        if (net_changed) {
            collect_path_occupancy_keys(segs, final_keys);
            replace_path_occupancy(self_keys, final_keys, occupied_set);
            changed_nets++;
        }
    }

    rebuild_occupied_from_paths(rc, routing_paths, occupied_set);
    std::cout << "[RunViaReduce] changed " << changed_runs
              << " runs in " << changed_nets
              << " nets, estimated saved vias " << saved_vias << std::endl;
    return saved_vias;
}

static bool path_has_sharp_angle(const std::vector<Seg>& segs) {
    struct IncidentDirection {
        ResultPoint vertex;
        ResultPoint other;
    };

    std::unordered_map<uint64_t, std::vector<IncidentDirection>> incidents;
    auto key_of = [](const ResultPoint& point) {
        return encode_coord(point.x, point.y, point.layer);
    };
    for (const Seg& seg : segs) {
        const ResultPoint& start = seg.first;
        const ResultPoint& end = seg.second;
        if (start.layer != end.layer ||
            (start.x == end.x && start.y == end.y)) {
            continue;
        }
        const int dx = (end.x > start.x) - (end.x < start.x);
        const int dy = (end.y > start.y) - (end.y < start.y);
        const int steps = std::max(
            std::abs(end.x - start.x), std::abs(end.y - start.y));
        ResultPoint previous = start;
        for (int step = 1; step <= steps; ++step) {
            ResultPoint current{
                start.x + dx * step, start.y + dy * step, start.layer};
            incidents[key_of(previous)].push_back(
                IncidentDirection{previous, current});
            incidents[key_of(current)].push_back(
                IncidentDirection{current, previous});
            previous = current;
        }
    }

    for (const auto& entry : incidents) {
        const std::vector<IncidentDirection>& arms = entry.second;
        for (size_t first = 0; first < arms.size(); ++first) {
            for (size_t second = first + 1; second < arms.size(); ++second) {
                const double angle = angle_between_pts(
                    arms[first].other, arms[first].vertex,
                    arms[second].other);
                if (angle > 0.0 && angle < 90.0) return true;
            }
        }
    }
    return false;
}

void smooth_sharp_angles(
    std::unordered_map<std::string, std::vector<Seg>>& routing_paths,
    OccGrid& occupied_set,
    int layer_count,
    int routing_threads
) {
    int fixed = 0, skipped = 0;
    std::unordered_set<std::string> candidate_nets;
    if (routing_threads != 1 && !routing_paths.empty()) {
        std::vector<const std::pair<const std::string, std::vector<Seg>>*>
            path_entries;
        path_entries.reserve(routing_paths.size());
        for (const auto& entry : routing_paths) path_entries.push_back(&entry);

        std::vector<uint8_t> has_sharp_angle(path_entries.size(), 0);
        const int screening_threads = routing_threads > 0
            ? std::min<int>(routing_threads, path_entries.size())
            : std::min<int>(omp_get_max_threads(), path_entries.size());
#pragma omp parallel for schedule(dynamic, 1) num_threads(screening_threads)
        for (size_t index = 0; index < path_entries.size(); ++index) {
            has_sharp_angle[index] = path_has_sharp_angle(
                path_entries[index]->second);
        }
        for (size_t index = 0; index < path_entries.size(); ++index) {
            if (has_sharp_angle[index]) {
                candidate_nets.insert(path_entries[index]->first);
            }
        }
        std::cout << "[SmoothScreen] candidates " << candidate_nets.size()
                  << " of " << routing_paths.size()
                  << ", threads " << screening_threads << std::endl;
    }
    // 多趟：每趟每个 net 修一处锐角，重复直到无可修（或达上限），覆盖一个 net 多个锐角
    skipped = 0;
    for (auto& kv : routing_paths) {
        if (routing_threads != 1 &&
            candidate_nets.find(kv.first) == candidate_nets.end()) {
            continue;
        }
        for (int pass = 0; pass < 30; ++pass) {
        int fixed_this_pass = 0;
        std::vector<Seg>& segs = kv.second;
        // Expand planar segments into unit incident edges, matching evaluator.py.
        struct IncidentArm {
            size_t segment_index;
            ResultPoint vertex;
            ResultPoint other;
        };
        std::unordered_map<uint64_t, std::vector<IncidentArm>> inc;
        auto key_of = [](const ResultPoint& p) { return encode_coord(p.x, p.y, p.layer); };
        std::unordered_set<uint64_t> self_points;
        self_points.reserve(segs.size() * 2);
        for (const Seg& seg : segs) {
            self_points.insert(key_of(seg.first));
            self_points.insert(key_of(seg.second));
        }
        for (size_t i = 0; i < segs.size(); ++i) {
            const ResultPoint& a = segs[i].first;
            const ResultPoint& b = segs[i].second;
            if (a.layer != b.layer || (a.x == b.x && a.y == b.y)) continue;
            int dx = (b.x > a.x) - (b.x < a.x);
            int dy = (b.y > a.y) - (b.y < a.y);
            int steps = std::max(std::abs(b.x - a.x), std::abs(b.y - a.y));
            ResultPoint previous = a;
            for (int step = 1; step <= steps; ++step) {
                ResultPoint current{a.x + dx * step, a.y + dy * step, a.layer};
                inc[key_of(previous)].push_back(IncidentArm{i, previous, current});
                inc[key_of(current)].push_back(IncidentArm{i, current, previous});
                previous = current;
            }
        }
        for (auto& vk : inc) {
            std::vector<IncidentArm>& arms = vk.second;
            if (arms.size() < 2) continue;
            ResultPoint B = arms[0].vertex;
            size_t pi = SIZE_MAX, pj = SIZE_MAX;
            ResultPoint Fi, Fj;
            for (size_t a = 0; a < arms.size() && pi == SIZE_MAX; ++a) {
                const ResultPoint& Fa = arms[a].other;
                for (size_t b = a + 1; b < arms.size(); ++b) {
                    const ResultPoint& Fb = arms[b].other;
                    double ang = angle_between_pts(Fa, B, Fb);
                    if (ang > 0.0 && ang < 90.0) {
                        pi = arms[a].segment_index;
                        pj = arms[b].segment_index;
                        Fi = Fa;
                        Fj = Fb;
                        break;
                    }
                }
            }
            if (pi == SIZE_MAX) continue;
            size_t i0 = pi, i1 = pj;
            ResultPoint F0 = Fi, F1 = Fj;
            {
            // 选较短臂抬升（按切比雪夫距离）
            int d0 = std::max(std::abs(F0.x - B.x), std::abs(F0.y - B.y));
            int d1 = std::max(std::abs(F1.x - B.x), std::abs(F1.y - B.y));
            size_t lift_candidates[2] = {i0, i1};
            ResultPoint endpoint_candidates[2] = {F0, F1};
            if (d1 < d0) {
                std::swap(lift_candidates[0], lift_candidates[1]);
                std::swap(endpoint_candidates[0], endpoint_candidates[1]);
            }
            int L = B.layer;
            auto find_layer = [&](const ResultPoint& endpoint) {
                for (int cand = L + 1; cand < layer_count; ++cand) {
                    if (via_dest_free(B.x, B.y, cand, occupied_set, self_points, segs) &&
                        via_dest_free(endpoint.x, endpoint.y, cand, occupied_set, self_points, segs)) {
                        return cand;
                    }
                }
                for (int cand = L - 1; cand >= 0; --cand) {
                    if (via_dest_free(B.x, B.y, cand, occupied_set, self_points, segs) &&
                        via_dest_free(endpoint.x, endpoint.y, cand, occupied_set, self_points, segs)) {
                        return cand;
                    }
                }
                return -1;
            };

            size_t lift = SIZE_MAX;
            ResultPoint F;
            int L1 = -1;
            for (int candidate = 0; candidate < 2; ++candidate) {
                int candidate_layer = find_layer(endpoint_candidates[candidate]);
                if (candidate_layer < 0) continue;
                lift = lift_candidates[candidate];
                F = endpoint_candidates[candidate];
                L1 = candidate_layer;
                break;
            }
            if (L1 < 0) { skipped++; continue; }
            // Replace only the selected unit edge; preserve both sides of a long segment.
            Seg original = segs[lift];
            int step_x = (original.second.x > original.first.x) -
                         (original.second.x < original.first.x);
            int step_y = (original.second.y > original.first.y) -
                         (original.second.y < original.first.y);
            bool forward = (F.x == B.x + step_x && F.y == B.y + step_y);
            std::vector<Seg> replacement_segments;
            auto add_planar = [&](const ResultPoint& from, const ResultPoint& to) {
                if (from.x != to.x || from.y != to.y || from.layer != to.layer) {
                    replacement_segments.emplace_back(from, to);
                }
            };
            if (forward) {
                add_planar(original.first, B);
                add_planar(F, original.second);
            } else {
                add_planar(original.first, F);
                add_planar(B, original.second);
            }
            replacement_segments.emplace_back(
                ResultPoint{B.x, B.y, L}, ResultPoint{B.x, B.y, L1});
            replacement_segments.emplace_back(
                ResultPoint{B.x, B.y, L1}, ResultPoint{F.x, F.y, L1});
            replacement_segments.emplace_back(
                ResultPoint{F.x, F.y, L1}, ResultPoint{F.x, F.y, L});
            segs[lift] = replacement_segments.front();
            for (size_t replacement_idx = 1;
                 replacement_idx < replacement_segments.size(); ++replacement_idx) {
                segs.push_back(replacement_segments[replacement_idx]);
            }
            // 占用新增点：仅 L1 层两端点 + 上层平移段对角中点（与评测 occupancy 口径一致）
            occupied_set.insert(encode_coord(B.x, B.y, L1));
            occupied_set.insert(encode_coord(F.x, F.y, L1));
            if (B.x != F.x && B.y != F.y) {
                float mx = (B.x + F.x) * 0.5f, my = (B.y + F.y) * 0.5f;
                occupied_set.insert(encode_coord(mx, my, L1));
            }
            fixed++;
            fixed_this_pass++;
            break;  // 该 net 本趟只改一处，避免下标失配；多锐角靠多趟
            }
        }
        if (fixed_this_pass == 0) break;  // 本趟无改动 → 收敛
        }
    }
    std::cout << "[Smooth] fixed " << fixed << " sharp angles, skipped " << skipped << std::endl;
}

void output_routing_result(
    const RoutingCase& rc,
    const std::unordered_map<std::string, std::vector<std::pair<ResultPoint, ResultPoint>>>& routing_paths,
    const std::string& output_filename
) {
    json result;
    for (const std::string& net_name : rc.net_name_list) {
        json net_segments = json::array();
        for (const auto& raw_segment : routing_paths.at(net_name)) {
            json segment;
            segment["start_grid_coordinate"] = coord_to_json(raw_segment.first.x, raw_segment.first.y, raw_segment.first.layer, rc.layer_count);
            segment["end_grid_coordinate"] = coord_to_json(raw_segment.second.x, raw_segment.second.y, raw_segment.second.layer, rc.layer_count);
            net_segments.push_back(segment);
        }
        result[net_name] = std::move(net_segments);
    }

    std::ofstream ofs(output_filename);
    if (ofs.is_open()) {
        ofs << std::setw(4) << result << std::endl;
        std::cout << "\nRouting result saved to: " << output_filename << std::endl;
    } else {
        std::cerr << "\nError: Failed to open output file " << output_filename << std::endl;
    }
}
