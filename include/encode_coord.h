#ifndef ENCODE_COORD_H
#define ENCODE_COORD_H

#include <cstdint>
#include <cmath>
#include <string>

// 整数点和对角线半整数中点共用同一个双倍坐标系：
//   整数 (x, y)       -> (2*x + 20000, 2*y + 20000) （偶数）
//   半整数 (x, y)   -> (round(2*x) + 20000, round(2*y) + 20000) （奇数）
// 这样不依赖原坐标小于 10000，也不需要从 key 值域猜测点类型。
// 每个双倍坐标和 layer 各占 20 bit。
inline uint64_t encode_coord(int x, int y, int layer) {
    const int OFFSET2 = 20000;
    const int ix = 2 * x + OFFSET2;
    const int iy = 2 * y + OFFSET2;
    uint64_t key = 0;
    key |= (static_cast<uint64_t>(ix) & 0xFFFFF) << 40;
    key |= (static_cast<uint64_t>(iy) & 0xFFFFF) << 20;
    key |= (static_cast<uint64_t>(layer) & 0xFFFFF);
    return key;
}

inline uint64_t encode_coord(float x, float y, int layer) {
    const int OFFSET2 = 20000;
    int ix = static_cast<int>(std::round(x * 2.0f)) + OFFSET2;
    int iy = static_cast<int>(std::round(y * 2.0f)) + OFFSET2;
    uint64_t key = 0;
    key |= (static_cast<uint64_t>(ix) & 0xFFFFF) << 40;
    key |= (static_cast<uint64_t>(iy) & 0xFFFFF) << 20;
    key |= (static_cast<uint64_t>(layer) & 0xFFFFF);
    return key;
}

#endif
