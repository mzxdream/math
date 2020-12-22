#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include <mzx/math/math_util.h>

using namespace mzx;

inline void TestMathUtil(int count)
{
    auto a1 = MathUtilF::Epsilon();
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