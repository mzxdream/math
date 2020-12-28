#ifndef __MZX_TEST_FIXED_INL_H__
#define __MZX_TEST_FIXED_INL_H__

#include <cstdint>
#include <cmath>
#include <random>

using Fixed64 = FixedNumber<int64_t, 32>;
using MathUtilF64 = MathUtil<Fixed64>;

static constexpr float FIXED_EPS = 0.000001f;
static constexpr float FIXED_FMIN = -100000.0f;
static constexpr float FIXED_FMAX = 100000.0f;
static constexpr int64_t FIXED_IMIN = -100000;
static constexpr int64_t FIXED_IMAX = 100000;

inline int64_t RandInt(int64_t mina, int64_t maxa)
{
    auto f = static_cast<float>(rand()) / RAND_MAX;
    auto t = mina + static_cast<int64_t>((maxa - mina) * f);
    return std::max(mina, std::min(t, maxa));
}

inline float RandFloat(float mina, float maxa)
{
    auto f = static_cast<float>(rand()) / RAND_MAX;
    auto t = mina + f * (maxa - mina);
    return std::max(mina, std::min(t, maxa));
}

#endif