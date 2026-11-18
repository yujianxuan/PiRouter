#pragma once

#include <algorithm>
#include <climits>
#include <cstdint>
#include <memory>
#include <vector>

// Per-search best-g/closed state stored in lazily allocated 64x64 pages.
// Pages are retained per thread and cleared logically with a generation stamp.
class AStarStateGrid {
public:
    static constexpr int BLOCK_SHIFT = 6;
    static constexpr int BLOCK_SIZE = 1 << BLOCK_SHIFT;
    static constexpr int BLOCK_MASK = BLOCK_SIZE - 1;

    AStarStateGrid() = default;
    AStarStateGrid(const AStarStateGrid&) = delete;
    AStarStateGrid& operator=(const AStarStateGrid&) = delete;

    void reset(int width, int height, int layers) {
        width = std::max(width, 1);
        height = std::max(height, 1);
        layers = std::max(layers, 1);

        const int required_width_blocks =
            (width + BLOCK_SIZE - 1) >> BLOCK_SHIFT;
        const int required_height_blocks =
            (height + BLOCK_SIZE - 1) >> BLOCK_SHIFT;

        if (layers_ != layers ||
            required_width_blocks > width_blocks_ ||
            required_height_blocks > height_blocks_) {
            width_blocks_ = std::max(width_blocks_, required_width_blocks);
            height_blocks_ = std::max(height_blocks_, required_height_blocks);
            layers_ = layers;
            table_.clear();
            table_.resize(
                static_cast<size_t>(width_blocks_) * height_blocks_);
        }

        width_ = width;
        height_ = height;
        ++generation_;
        if (generation_ == 0) {
            clear_stamps();
            generation_ = 1;
        }
    }

    void set_best(int x, int y, int layer, int g) {
        Cell* cell = get_cell(x, y, layer, true);
        if (cell == nullptr) return;
        cell->generation = generation_;
        cell->best_g = g;
        cell->closed = false;
    }

    int best_g(int x, int y, int layer) const {
        const Cell* cell = get_cell(x, y, layer);
        if (cell == nullptr || cell->generation != generation_) return INT_MAX;
        return cell->best_g;
    }

    bool is_closed(int x, int y, int layer) const {
        const Cell* cell = get_cell(x, y, layer);
        return cell != nullptr && cell->generation == generation_ && cell->closed;
    }

    void close(int x, int y, int layer) {
        Cell* cell = get_cell(x, y, layer, true);
        if (cell == nullptr) return;
        if (cell->generation != generation_) {
            cell->generation = generation_;
            cell->best_g = INT_MAX;
        }
        cell->closed = true;
    }

    bool improve(int x, int y, int layer, int g) {
        Cell* cell = get_cell(x, y, layer, true);
        if (cell == nullptr) return false;
        if (cell->generation != generation_) {
            cell->generation = generation_;
            cell->best_g = g;
            cell->closed = false;
            return true;
        }
        if (g >= cell->best_g) return false;
        cell->best_g = g;
        return true;
    }

private:
    struct Cell {
        int best_g = INT_MAX;
        uint32_t generation = 0;
        bool closed = false;
    };

    struct Block {
        explicit Block(int layers)
            : cells(static_cast<size_t>(BLOCK_SIZE) * BLOCK_SIZE * layers) {}
        std::vector<Cell> cells;
    };

    Cell* get_cell(int x, int y, int layer, bool create) {
        if (x < 0 || x >= width_ || y < 0 || y >= height_ ||
            layer < 0 || layer >= layers_) {
            return nullptr;
        }
        const int block_index =
            (y >> BLOCK_SHIFT) * width_blocks_ + (x >> BLOCK_SHIFT);
        std::unique_ptr<Block>& block = table_[block_index];
        if (!block) {
            if (!create) return nullptr;
            block = std::make_unique<Block>(layers_);
        }
        const size_t cell_index =
            (static_cast<size_t>(x & BLOCK_MASK) * BLOCK_SIZE +
             (y & BLOCK_MASK)) * layers_ + layer;
        return &block->cells[cell_index];
    }

    const Cell* get_cell(int x, int y, int layer) const {
        if (x < 0 || x >= width_ || y < 0 || y >= height_ ||
            layer < 0 || layer >= layers_) {
            return nullptr;
        }
        const int block_index =
            (y >> BLOCK_SHIFT) * width_blocks_ + (x >> BLOCK_SHIFT);
        const std::unique_ptr<Block>& block = table_[block_index];
        if (!block) return nullptr;
        const size_t cell_index =
            (static_cast<size_t>(x & BLOCK_MASK) * BLOCK_SIZE +
             (y & BLOCK_MASK)) * layers_ + layer;
        return &block->cells[cell_index];
    }

    void clear_stamps() {
        for (const auto& block : table_) {
            if (!block) continue;
            for (Cell& cell : block->cells) cell.generation = 0;
        }
    }

    int width_ = 0;
    int height_ = 0;
    int layers_ = 0;
    int width_blocks_ = 0;
    int height_blocks_ = 0;
    uint32_t generation_ = 0;
    std::vector<std::unique_ptr<Block>> table_;
};
