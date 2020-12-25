#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include <mzx/math/math_util.h>

using namespace mzx;

template <typename T>
inline void CheckOverflow(T x, T y)
{
    if (IsAddOverflow(x, y))
    {
        std::cout << "x:" << x << " y:" << y << " add overflow" << std::endl;
    }
    else
    {
        std::cout << "x:" << x << " y:" << y << " add ok" << std::endl;
    }
    if (IsSubOverflow(x, y))
    {
        std::cout << "x:" << x << " y:" << y << " sub overflow" << std::endl;
    }
    else
    {
        std::cout << "x:" << x << " y:" << y << " sub ok" << std::endl;
    }
    if (IsMulOverflow(x, y))
    {
        std::cout << "x:" << x << " y:" << y << " mul overflow" << std::endl;
    }
    else
    {
        std::cout << "x:" << x << " y:" << y << " mul ok" << std::endl;
    }
    if (IsDivOverflow(x, y))
    {
        std::cout << "x:" << x << " y:" << y << " div overflow" << std::endl;
    }
    else
    {
        std::cout << "x:" << x << " y:" << y << " div ok" << std::endl;
    }
}

inline void TestMathUtil(int count)
{
    int64_t arr1[] = {-3, -2, -1, 0, 1, 2, 3, INT64_MIN / 3, INT64_MIN / 2, INT64_MIN + 1, INT64_MIN, INT64_MAX / 3, INT64_MAX / 2, INT64_MAX - 1, INT64_MAX};
    for (auto x : arr1)
    {
        for (auto y : arr1)
        {
            CheckOverflow(x, y);
        }
    }
    uint64_t arr2[] = {0, 1, 2, 3, UINT64_MAX / 3, UINT64_MAX / 2, UINT64_MAX};
    for (auto x : arr2)
    {
        for (auto y : arr2)
        {
            CheckOverflow(x, y);
        }
    }

    auto a2 = MathUtilF::PI();
    auto a3 = MathUtilF::HalfPI();
    auto a4 = MathUtilF::TwoPI();
    auto a5 = MathUtilF::Sqrt2();
    auto a6 = MathUtilF::HalfSqrt2();
    auto a7 = MathUtilF::TwoSqrt2();

    auto b1 = MathUtilF::CastFrom(123);
    auto b2 = MathUtilF::CastFrom(1, 2);
    auto b3 = MathUtilF::CompareApproximately(1, 1.1);
    auto b4 = MathUtilF::Abs(-1.1);
    auto b5 = MathUtilF::Min(1, 2);
    auto b6 = MathUtilF::Max(1, 2);
    auto b7 = MathUtilF::Clamp(1, 2, 3);
    auto b8 = MathUtilF::Sqrt(2);
    auto b9 = MathUtilF::Cos(12);
    auto b10 = MathUtilF::Sin(12);
    auto b11 = MathUtilF::Acos(0.5);
    auto b12 = MathUtilF::Asin(0.5);
    auto b13 = MathUtilF::Atan(1);
    auto b14 = MathUtilF::Atan2(2, 3);
    auto b15 = MathUtilF::Rad2Deg(12);
    auto b16 = MathUtilF::Deg2Rad(12);
}