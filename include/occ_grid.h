#ifndef OCC_GRID_H
#define OCC_GRID_H

#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>

// ============================================================================
// OccGrid: 用分页稀疏网格替代 std::unordered_set<uint64_t> 做占用集合。
//   - 64x64 分块（block），block 触碰才分配，内存稀疏。
//   - 数组下标 O(1) 访问，连续内存 cache 友好，无哈希/rehash。
//   - encode_coord 已将整数点与半整数中点统一到双倍分辨率：
//       整数点 (x,y)          -> 网格下标 (2*x+20000, 2*y+20000)
//       半整数中点 (x+0.5...)  -> 网格下标 (奇数, 奇数)
//     两者落在不同奇偶下标，天然不冲突，无需单独容器。
//   - 接口 count(key)/insert(key) 与原 encode_coord 键语义完全一致，
//     调用点只需把 `.find(k)!=.end()` 换成 `.count(k)`。
// ============================================================================

class OccGrid {
public:
    static const int BLOCK = 64;
    static const int SHIFT = 6;     // 2^6 = 64
    static const int MASK = 0x3F;   // 63

private:
    struct Block {
        // gx,gy 取 block 内偏移；layer 维度运行时大小
        std::vector<uint8_t> data;  // [BLOCK*BLOCK*layers]
        explicit Block(int layers) : data((size_t)BLOCK * BLOCK * layers, 0) {}
    };

    int layers_ = 0;
    int gw_ = 0, gh_ = 0;            // 双倍分辨率网格尺寸
    int wblocks_ = 0, hblocks_ = 0;
    std::vector<Block*> table_;
    std::vector<int> allocated_;     // 已分配 block 的 table 下标，供快速遍历/合并

    inline int block_idx(int gx, int gy) const {
        return (gy >> SHIFT) * wblocks_ + (gx >> SHIFT);
    }

public:
    OccGrid() = default;

    // 双倍分辨率：原始坐标范围用 OFFSET_INT 偏移后约 [0, 2*OFFSET_INT]，×2 → 网格边
    void init(int max_coord_with_offset, int layers) {
        layers_ = layers;
        gw_ = gh_ = 2 * (max_coord_with_offset + 2);  // +2 余量给奇数中点
        wblocks_ = (gw_ + BLOCK - 1) / BLOCK;
        hblocks_ = (gh_ + BLOCK - 1) / BLOCK;
        table_.assign((size_t)wblocks_ * hblocks_, nullptr);
        allocated_.clear();
    }

    ~OccGrid() { clear_blocks(); }

    void clear_blocks() {
        for (int idx : allocated_) { delete table_[idx]; table_[idx] = nullptr; }
        allocated_.clear();
    }

    // --- 从 encode_coord 的 64 位 key 解出双倍网格下标 (gx,gy,layer) ---
    // key 中 x/y 字段已是统一双倍坐标，不再按值域猜测整数/半整数类型。
    inline void decode(uint64_t key, int& gx, int& gy, int& layer) const {
        gx = (int)((key >> 40) & 0xFFFFF);
        gy = (int)((key >> 20) & 0xFFFFF);
        layer = (int)(key & 0xFFFFF);
    }

    inline bool count(uint64_t key) const {
        if (table_.empty()) return false;  // 未初始化保护
        int gx, gy, l;
        decode(key, gx, gy, l);
        if (gx < 0 || gx >= gw_ || gy < 0 || gy >= gh_ || l < 0 || l >= layers_) return false;
        int idx = block_idx(gx, gy);
        if (idx < 0 || idx >= (int)table_.size()) return false;
        Block* b = table_[idx];
        if (!b) return false;
        size_t off = ((size_t)(gx & MASK) * BLOCK + (gy & MASK)) * layers_ + l;
        return b->data[off] != 0;
    }

    // Check the same lattice points used by the evaluator for one planar
    // queen-line segment. Integer horizontal/vertical routes visit even
    // doubled coordinates; diagonals additionally visit odd midpoint cells.
    // Empty sparse blocks are skipped as a unit, which is important for the
    // long pattern candidates in large interposer designs.
    inline bool planar_segment_clear(
        int x0, int y0, int x1, int y1, int layer,
        bool exclude_final = false) const {

        const int dx = (x1 > x0) - (x1 < x0);
        const int dy = (y1 > y0) - (y1 < y0);
        const int span_x = std::abs(x1 - x0);
        const int span_y = std::abs(y1 - y0);
        if (dx != 0 && dy != 0 && span_x != span_y) return false;
        const bool diagonal = dx != 0 && dy != 0;
        const int total_steps = diagonal
            ? 2 * span_x : std::max(span_x, span_y);
        if (total_steps == 0) return true;

        constexpr int OFFSET2 = 20000;
        int gx = 2 * x0 + OFFSET2;
        int gy = 2 * y0 + OFFSET2;
        const int step_x = diagonal ? dx : 2 * dx;
        const int step_y = diagonal ? dy : 2 * dy;
        int step = 1;
        gx += step_x;
        gy += step_y;
        while (step <= total_steps) {
            if (!(exclude_final && step == total_steps)) {
                if (gx < 0 || gx >= gw_ || gy < 0 || gy >= gh_ ||
                    layer < 0 || layer >= layers_) {
                    return false;
                }
                const int idx = block_idx(gx, gy);
                Block* block = table_[idx];
                if (block != nullptr) {
                    const size_t offset =
                        ((size_t)(gx & MASK) * BLOCK + (gy & MASK)) *
                        layers_ + layer;
                    if (block->data[offset] != 0) return false;
                } else {
                    int jump_x = total_steps - step + 1;
                    int jump_y = jump_x;
                    if (step_x > 0) {
                        jump_x = (((gx >> SHIFT) + 1) * BLOCK - gx +
                                  step_x - 1) / step_x;
                    } else if (step_x < 0) {
                        jump_x = (gx - ((gx >> SHIFT) * BLOCK) +
                                  (-step_x) - 1) / (-step_x) + 1;
                    }
                    if (step_y > 0) {
                        jump_y = (((gy >> SHIFT) + 1) * BLOCK - gy +
                                  step_y - 1) / step_y;
                    } else if (step_y < 0) {
                        jump_y = (gy - ((gy >> SHIFT) * BLOCK) +
                                  (-step_y) - 1) / (-step_y) + 1;
                    }
                    const int jump = std::max(
                        1, std::min({jump_x, jump_y,
                                     total_steps - step + 1}));
                    step += jump;
                    gx += jump * step_x;
                    gy += jump * step_y;
                    continue;
                }
            }
            ++step;
            gx += step_x;
            gy += step_y;
        }
        return true;
    }

    inline void insert(uint64_t key) {
        if (table_.empty()) return;  // 未初始化保护
        int gx, gy, l;
        decode(key, gx, gy, l);
        if (gx < 0 || gx >= gw_ || gy < 0 || gy >= gh_ || l < 0 || l >= layers_) return;
        int idx = block_idx(gx, gy);
        if (idx < 0 || idx >= (int)table_.size()) return;
        Block* b = table_[idx];
        if (!b) { b = new Block(layers_); table_[idx] = b; allocated_.push_back(idx); }
        size_t off = ((size_t)(gx & MASK) * BLOCK + (gy & MASK)) * layers_ + l;
        b->data[off] = 1;
    }

    inline void erase(uint64_t key) {
        if (table_.empty()) return;
        int gx, gy, l;
        decode(key, gx, gy, l);
        if (gx < 0 || gx >= gw_ || gy < 0 || gy >= gh_ || l < 0 || l >= layers_) return;
        int idx = block_idx(gx, gy);
        if (idx < 0 || idx >= (int)table_.size()) return;
        Block* b = table_[idx];
        if (!b) return;
        size_t off = ((size_t)(gx & MASK) * BLOCK + (gy & MASK)) * layers_ + l;
        b->data[off] = 0;
    }

    // 多线程合并：把 other 的所有已占格点并入本网格
    void merge_from(const OccGrid& other) {
        for (int idx : other.allocated_) {
            Block* ob = other.table_[idx];
            if (!ob) continue;
            Block* b = table_[idx];
            if (!b) { b = new Block(layers_); table_[idx] = b; allocated_.push_back(idx); }
            const std::vector<uint8_t>& od = ob->data;
            std::vector<uint8_t>& d = b->data;
            for (size_t i = 0; i < od.size(); ++i) if (od[i]) d[i] = 1;
        }
    }

    // --- STL 兼容接口：让 `occ.find(k) != occ.end()` 零改动可用 ---
    // find() 返回一个布尔代理；end() 返回 false 代理；二者用 != 比较即"是否存在"。
    struct Iter {
        bool present;
        bool operator!=(const Iter& o) const { return present != o.present; }
        bool operator==(const Iter& o) const { return present == o.present; }
    };
    inline Iter find(uint64_t key) const { return Iter{ count(key) }; }
    inline Iter end() const { return Iter{ false }; }
};

#endif
