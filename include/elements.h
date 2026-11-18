#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
#include "encode_coord.h"
#include "occ_grid.h"

struct ResultPoint {
    int x;
    int y;
    int layer;
};

typedef std::pair<ResultPoint, ResultPoint> Seg;
typedef std::unordered_map<std::string, std::vector<Seg>> DRResultDict, ResultDict;

struct ExitTerminal {
    std::string bump_name;
    int grid_coord_x;
    int grid_coord_y;
    int layer;

    bool operator==(const ExitTerminal& other) const {
        return grid_coord_x == other.grid_coord_x
               && grid_coord_y == other.grid_coord_y
               && layer == other.layer;
    }
};

struct Bump {
    std::string net_name;
    int bump_type;
    int grid_coord_x;
    int grid_coord_y;
    int layer;
    ExitTerminal exit_terminal;
    size_t pair_index = std::numeric_limits<size_t>::max();

    bool operator==(const Bump& other) const {
        return net_name == other.net_name
               && bump_type == other.bump_type
               && grid_coord_x == other.grid_coord_x
               && grid_coord_y == other.grid_coord_y
               && layer == other.layer;
    }
};

struct Net {
    int net_type;
    int bump_count;
    std::vector<std::string> bump_names;
};

struct RoutingCase {
    std::string case_name;
    std::string output_dir;
    int grid_length;
    int grid_max_width;
    int grid_max_height;
    int layer_count;
    std::unordered_map<std::string, Bump> bump_dict;
    int net_count;
    std::vector<std::string> net_name_list;
    std::unordered_map<std::string, Net> net_dict;
    OccGrid occupied_set;
    std::unordered_set<uint64_t> occupied_set_2d;
};
