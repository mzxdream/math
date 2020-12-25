#include <iostream>
#include <mzx/math/fixed_math_util.h>
#include "test_fixed_inl.h"

using namespace mzx;
using Fixed64 = FixedNumber<int64_t, 32>;
using MathUtilF64 = MathUtil<Fixed64>;

static void TestAbs(int count)
{
    auto diff_i = 0;
    auto diff_a = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = a >= 0 ? a : -a;
        auto t2 = MathUtilF64::Abs(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_i = i;
            diff_max = diff;
            diff_a = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestAbs i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestAbs count:" << count << " diff i:" << diff_i << " diff max:" << diff_max << " a:" << diff_a << " diff count:" << diff_count << std::endl;
}

static void TestMin(int count)
{
    auto diff_i = 0;
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto b = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = std::min(a, b);
        auto t2 = MathUtilF64::Min(Fixed64::FromFloat(a), Fixed64::FromFloat(b));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_i = i;
            diff_max = diff;
            diff_a = a;
            diff_b = b;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestMin i:" << i << " a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestMin count:" << count << " diff i:" << diff_i << " diff max:" << diff_max << " a:" << diff_a << " b:" << diff_b << " diff count:" << diff_count << std::endl;
}

static void TestMax(int count)
{
    auto diff_i = 0;
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto b = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = std::max(a, b);
        auto t2 = MathUtilF64::Max(Fixed64::FromFloat(a), Fixed64::FromFloat(b));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_i = i;
            diff_max = diff;
            diff_a = a;
            diff_b = b;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestMax i:" << i << " a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestMax count:" << count << " diff i:" << diff_i << " diff max:" << diff_max << " a:" << diff_a << " b:" << diff_b << " diff count:" << diff_count << std::endl;
}

static void TestClamp(int count)
{
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_c = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto b = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto c = RandFloat(FIXED_FMIN, FIXED_FMAX);
        if (b > c)
        {
            std::swap(b, c);
        }
        auto t1 = std::min(std::max(a, b), c);
        auto t2 = MathUtilF64::Clamp(Fixed64::FromFloat(a), Fixed64::FromFloat(b), Fixed64::FromFloat(c));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_a = a;
            diff_b = b;
            diff_c = c;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestClamp i:" << i << " a:" << t1
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestClamp count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " b:" << diff_b << " c:" << diff_c << " diff count:" << diff_count << std::endl;
}

static void TestSqrt(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(0.0001, FIXED_FMAX);
        auto t1 = sqrt(a);
        auto t2 = MathUtilF64::Sqrt(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestSqrt i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestSqrt count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestPow(int count)
{
    auto diff_i = 0;
    auto diff_a = 0.0;
    auto diff_n = 0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto n = static_cast<int>(RandInt(-5, 5));
        auto t1 = pow(a, n);
        if (abs(t1) >= FIXED_FMAX - 1)
        {
            --i;
            continue;
        }
        auto t2 = MathUtilF64::Pow(Fixed64::FromFloat(a), n);
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_i = i;
            diff_max = diff;
            diff_a = a;
            diff_n = n;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestPow i:" << i << " a:" << a << " n:" << n
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestPow count:" << count << " diff i:" << diff_i << " diff max:" << diff_max << " a:" << diff_a << " n:" << diff_n << " diff count:" << diff_count << std::endl;
}

static void TestCosDeg(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX) / 2;
        auto t1 = cos(a * 3.14159265358979323846 / 180);
        auto t2 = MathUtilF64::CosDeg(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestCosDeg i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestCosDeg count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestSinDeg(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX) / 2;
        auto t1 = sin(a * 3.14159265358979323846 / 180);
        auto t2 = MathUtilF64::SinDeg(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestSinDeg i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestSinDeg count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestCos(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = cos(a);
        auto t2 = MathUtilF64::Cos(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestCos i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << t1 - t2.ToFloat() << std::endl;
        }
    }
    std::cout << "TestCos count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestSin(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = sin(a);
        auto t2 = MathUtilF64::Sin(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestSin i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestSin count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestAcos(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(-1.0, 1.0);
        auto t1 = acos(a);
        auto t2 = MathUtilF64::Acos(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestAcos i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << t1 - t2.ToFloat() << std::endl;
        }
    }
    std::cout << "TestAcos count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestAsin(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(-1.0, 1.0);
        auto t1 = asin(a);
        auto t2 = MathUtilF64::Asin(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestAsin i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << t1 - t2.ToFloat() << std::endl;
        }
    }
    std::cout << "TestAsin count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestAtan(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = atan(a);
        auto t2 = MathUtilF64::Atan(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestAtan i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << t1 - t2.ToFloat() << std::endl;
        }
    }
    std::cout << "TestAtan count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestAtan2(int count)
{
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto b = RandFloat(FIXED_FMIN, FIXED_FMAX);
        if (abs(a) < 0.000001 && abs(b) < 0.000001)
        {
            --i;
            continue;
        }
        auto t1 = atan2(b, a);
        auto t2 = MathUtilF64::Atan2(Fixed64::FromFloat(b), Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_a = a;
            diff_b = b;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestAtan2 i:" << i << " a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << t1 - t2.ToFloat() << std::endl;
        }
    }
    std::cout << "TestAtan2 count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " b:" << diff_b << " diff count:" << diff_count << std::endl;
}

static void TestRad2Deg(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = a * 180 / 3.14159265358979323846;
        int64_t k = ((int64_t)(t1 / 360));
        int64_t k2 = k * 360;
        t1 -= ((int64_t)(t1 / 360)) * 360;
        while (t1 < 0)
        {
            t1 += 360;
        }
        auto t2 = MathUtilF64::Rad2Deg(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestRad2Deg i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestRad2Deg count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

static void TestDeg2Rad(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = a - ((int64_t)a / 360) * 360;
        if (t1 < 0)
        {
            t1 += 360;
        }
        t1 = t1 * 3.14159265358979323846 / 180;
        auto t2 = MathUtilF64::Deg2Rad(Fixed64::FromFloat(a));
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = a;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestDeg2Rad i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestDeg2Rad count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

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

    TestAbs(count);
    TestMin(count);
    TestMax(count);
    TestClamp(count);
    TestSqrt(count);
    TestPow(count);
    TestCosDeg(count);
    TestSinDeg(count);
    TestCos(count);
    TestSin(count);
    TestAcos(count);
    TestAsin(count);
    TestAtan(count);
    TestAtan2(count);
    TestRad2Deg(count);
    TestDeg2Rad(count);
}