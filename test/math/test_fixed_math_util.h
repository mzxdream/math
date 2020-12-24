#include <iostream>
#include <mzx/math/fixed_math_util.h>

using namespace mzx;
using Fixed64 = FixedNumber<int64_t, 32>;
using MathUtilF64 = MathUtil<Fixed64>;

void TestFixedMathUtil(int count)
{
    auto a1 = MathUtilF64::PI();
    auto a2 = MathUtilF64::HalfPI();
    auto a3 = MathUtilF64::TwoPI();
    auto a4 = MathUtilF64::Sqrt2();
    auto a5 = MathUtilF64::HalfSqrt2();
    auto a6 = MathUtilF64::TwoSqrt2();
    auto a7 = MathUtilF64::CastFrom(1);
    auto a8 = MathUtilF64::CastFrom(1, 2);
    auto a9 = MathUtilF64::CompareApproximately(Fixed64::FromFloat(1.1f), Fixed64::FromFloat(1.000001f));
    auto a10 = MathUtilF64::Abs(Fixed64::FromInt(-3));
    auto a11 = MathUtilF64::Min(Fixed64::FromInt(-1), Fixed64::FromInt(-2));
    auto a12 = MathUtilF64::Max(Fixed64::FromInt(-1), Fixed64::FromInt(-2));
    auto a13 = MathUtilF64::Clamp(Fixed64::FromInt(1), Fixed64::FromInt(2), Fixed64::FromFloat(3));
    auto a14 = MathUtilF64::Lerp(Fixed64::FromInt(-1), Fixed64::FromInt(-2), Fixed64::FromFloat(0.5f));
    auto a15 = MathUtilF64::Sqrt(Fixed64::FromInt(3));
    auto a16 = MathUtilF64::Pow(Fixed64::FromInt(3), -3);
    auto a17 = MathUtilF64::CosDeg(Fixed64::FromInt(180));
    auto a18 = MathUtilF64::SinDeg(Fixed64::FromInt(180));
    auto a19 = MathUtilF64::Cos(Fixed64::FromInt(1));
    auto a20 = MathUtilF64::Sin(Fixed64::FromInt(1));
    auto a21 = MathUtilF64::Acos(Fixed64::FromInt(1));
    auto a22 = MathUtilF64::Asin(Fixed64::FromInt(1));
    auto a23 = MathUtilF64::Atan(Fixed64::FromInt(1));
    auto a24 = MathUtilF64::Atan2(Fixed64::FromInt(2), Fixed64::FromInt(3));
    auto a26 = MathUtilF64::Rad2Deg(Fixed64::FromFloat(3.14f));
    auto a27 = MathUtilF64::Deg2Rad(Fixed64::FromFloat(90.0f));
}