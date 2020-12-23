#include <iostream>
#include <mzx/math/fixed_number.h>

using namespace mzx;

static const auto FIXED_EPS = 0.000001f;
static const auto FIXED_FMIN = -100000;
static const auto FIXED_FMAX = 100000;
static const auto FIXED_INT_MIN = -100000.0f;
static const auto FIXED_INT_MAX = 100000.0f;

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

void TestNaNInf(int count)
{
    float arr[] = {-1.0f, 1.0f, 0.0f, -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()};
    auto len = sizeof(arr) / sizeof(arr[0]);
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            auto k1 = arr[i];
            auto k2 = arr[j];

            auto p1 = k1 + k2;
            auto q1 = (Fixed64::FromFloat(k1) + Fixed64::FromFloat(k2)).ToFloat();
            auto p2 = k1 - k2;
            auto q2 = (Fixed64::FromFloat(k1) - Fixed64::FromFloat(k2)).ToFloat();
            auto p3 = k1 * k2;
            auto q3 = (Fixed64::FromFloat(k1) * Fixed64::FromFloat(k2)).ToFloat();
            auto p4 = k1 / k2;
            auto q4 = (Fixed64::FromFloat(k1) / Fixed64::FromFloat(k2)).ToFloat();
            auto p5 = -k1;
            auto q5 = (-Fixed64::FromFloat(k1)).ToFloat();
            auto p6 = k1 == k2;
            auto q6 = Fixed64::FromFloat(k1) == Fixed64::FromFloat(k2);
            auto p7 = k1 != k2;
            auto q7 = Fixed64::FromFloat(k1) != Fixed64::FromFloat(k2);
            auto p8 = k1 < k2;
            auto q8 = Fixed64::FromFloat(k1) < Fixed64::FromFloat(k2);
            auto p9 = k1 <= k2;
            auto q9 = Fixed64::FromFloat(k1) <= Fixed64::FromFloat(k2);
            auto p10 = k1 > k2;
            auto q10 = Fixed64::FromFloat(k1) > Fixed64::FromFloat(k2);
            auto p11 = k1 >= k2;
            auto q11 = Fixed64::FromFloat(k1) >= Fixed64::FromFloat(k2);

            std::cout << " i:" << i << " j:" << j << " k1:" << k1 << " k2:" << k2 << std::endl;
            std::cout << " p1:" << p1 << " p2:" << p2 << " p3:" << p3 << " p4:" << p4 << " p5:" << p5 << " p6:" << p6 << " p7:" << p7 << " p8:" << p8 << " p9:" << p9 << " p10:" << p10 << " p11:" << p11 << std::endl;
            std::cout << " q1:" << q1 << " q2:" << q2 << " q3:" << q3 << " q4:" << q4 << " q5:" << q5 << " q6:" << q6 << " q7:" << q7 << " q8:" << q8 << " q9:" << q9 << " q10:" << q10 << " q11:" << q11 << std::endl;
        }
    }
}

void TestFromInt(int count)
{
    int64_t diff_num = 0;
    auto diff_max = 0.0f;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto t1 = RandInt(FIXED_INT_MIN, FIXED_INT_MAX);
        auto t2 = Fixed64::FromInt(t1);
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestFromInt i:" << i << " a:" << t1
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestFromInt count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

void TestFromFloat(int count)
{
    auto diff_num = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto t1 = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t2 = Fixed64::FromFloat(t1);
        auto diff = abs(t1 - t2.ToFloat());
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestFromFloat i:" << i << " a:" << t1
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestFromFloat count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

void TestFromFraction(int count)
{
    int64_t diff_a = 0;
    int64_t diff_b = 0;
    auto diff_max = 0.0f;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandInt(FIXED_INT_MIN, FIXED_INT_MAX);
        auto b = RandInt(FIXED_INT_MIN, FIXED_INT_MAX);
        if (b == 0)
        {
            --i;
            continue;
        }
        auto t1 = static_cast<float>(a) / b;
        auto t2 = Fixed64::FromFraction(a, b);
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
            std::cout << "TestFromFraction i:" << i << " a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestFromFraction count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " b:" << diff_b << " diff count:" << diff_count << std::endl;
}

void TestToInt(int count)
{
    int64_t diff_num = 0;
    int64_t diff_max = 0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = static_cast<int64_t>(a);
        auto p = Fixed64::FromFloat(a);
        auto t2 = p.ToInt();
        auto diff = abs(t1 - t2);
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestToInt i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2 << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestToInt count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

void TestToFloorInt(int count)
{
    int64_t diff_num = 0;
    int64_t diff_max = 0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = static_cast<int64_t>(floor(a));
        auto p = Fixed64::FromFloat(a);
        auto t2 = p.ToFloorInt();
        auto diff = abs(t1 - t2);
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestToFloor i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2 << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestToFloor count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

void TestToCeilInt(int count)
{
    int64_t diff_num = 0;
    int64_t diff_max = 0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = static_cast<int64_t>(ceil(a));
        auto p = Fixed64::FromFloat(a);
        auto t2 = p.ToCeilInt();
        auto diff = abs(t1 - t2);
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestToCeil i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2 << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestToCeil count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

void TestToRoundInt(int count)
{
    int64_t diff_num = 0;
    int64_t diff_max = 0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = static_cast<int64_t>(round(a));
        auto p = Fixed64::FromFloat(a);
        auto t2 = p.ToRoundInt();
        auto diff = abs(t1 - t2);
        if (diff_max < diff)
        {
            diff_max = diff;
            diff_num = t1;
        }
        if (diff >= FIXED_EPS)
        {
            ++diff_count;
            std::cout << "TestToRound i:" << i << " a:" << a
                      << " t1:" << t1 << " t2:" << t2 << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestToRound count:" << count << " diff max:" << diff_max << " diff num:" << diff_num << " diff count:" << diff_count << std::endl;
}

void TestAdd(int count)
{
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto b = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = a + b;
        if (abs(t1) >= FIXED_FMAX - 1)
        {
            --i;
            continue;
        }
        auto t2 = Fixed64::FromFloat(a) + Fixed64::FromFloat(b);
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
            std::cout << "TestAdd a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestAdd count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " diff b:" << diff_b << " diff count:" << diff_count << std::endl;
}

void TestSub(int count)
{
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto b = RandFloat(FIXED_FMIN, FIXED_FMAX);
        auto t1 = a - b;
        if (abs(t1) >= FIXED_FMAX - 1)
        {
            --i;
            continue;
        }
        auto t2 = Fixed64::FromFloat(a) - Fixed64::FromFloat(b);
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
            std::cout << "TestSub a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestSub count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " diff b:" << diff_b << " diff count:" << diff_count << std::endl;
}

void TestMul(int count)
{
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN + 1, FIXED_FMAX - 1);
        auto b = RandFloat(FIXED_FMIN + 1, FIXED_FMAX - 1);
        a /= b;
        //a = -0.275013;
        //b = -5.48514e+11;
        //a = -2.7265845873233010;
        //b = -183967787690.66412;
        auto t1 = a * b;
        if (abs(t1) >= FIXED_FMAX - 1)
        {
            --i;
            continue;
        }
        auto k1 = Fixed64::FromFloat(a);
        auto k2 = Fixed64::FromFloat(b);
        t1 = k1.ToFloat() * k2.ToFloat();
        auto t2 = Fixed64::FromFloat(a) * Fixed64::FromFloat(b);
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
            std::cout << "TestMul a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestMul count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " diff b:" << diff_b << " diff count:" << diff_count << std::endl;
}

void TestDiv(int count)
{
    auto diff_a = 0.0;
    auto diff_b = 0.0;
    auto diff_max = 0.0;
    int diff_count = 0;
    for (int i = 0; i < count; i++)
    {
        auto a = RandFloat(FIXED_FMIN + 1, FIXED_FMAX - 1);
        auto b = RandFloat(FIXED_FMIN + 1, FIXED_FMAX - 1);
        b = a / b;
        //a = 1.67777e+07;
        //b = 3.63729e-05;
        //a = -157492532882.10339;
        //b = 0.28647724845118566;
        auto t1 = a / b;
        if (abs(t1) >= FIXED_FMAX / 2 - 1)
        {
            --i;
            continue;
        }
        auto k1 = Fixed64::FromFloat(a);
        auto k2 = Fixed64::FromFloat(b);
        t1 = k1.ToFloat() / k2.ToFloat();
        auto t2 = Fixed64::FromFloat(a) / Fixed64::FromFloat(b);
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
            std::cout << "TestDiv a:" << a << " b:" << b
                      << " t1:" << t1 << " t2:" << t2.ToFloat() << " diff:" << diff << std::endl;
        }
    }
    std::cout << "TestDiv count:" << count << " diff max:" << diff_max << " diff a:" << diff_a << " diff b:" << diff_b << " diff count:" << diff_count << std::endl;
}

inline void TestFixedNumber(int count)
{
    Fixed64 a(22222);
    Fixed64 b(234);
    a.Set(234);
    auto a1 = a.Get();
    auto a2 = a.ToFloat();
    auto a3 = a.ToInt();
    auto a4 = a.ToFloorInt();
    auto a5 = a.ToCeilInt();
    auto a6 = a.ToRoundInt();
    auto a7 = a.ToFloor();
    auto a8 = a.ToCeil();
    auto a9 = a.ToRound();
    auto a10 = a.IsNAN();
    auto a11 = a.IsPosInfinite();
    auto a12 = a.IsNegInfinite();
    auto a13 = a.IsFinite();
    auto a14 = a + b;
    auto a15 = a - b;
    auto a16 = -a;
    auto a17 = a * b;
    auto a18 = a / b;
    a += b;
    a -= b;
    a *= b;
    a /= b;
    auto a19 = a == b;
    auto a20 = a != b;
    auto a21 = a < b;
    auto a22 = a <= b;
    auto a23 = a > b;
    auto a24 = a >= b;
    auto a25 = Fixed64::FromInt(123);
    auto a26 = Fixed64::FromFloat(123.456);
    auto a27 = Fixed64::FromFraction(1, 2);
    auto a28 = Fixed64::Nan();
    auto a29 = Fixed64::Infinity();
    auto a30 = Fixed64::Epsilon();
    auto a31 = Fixed64::Min();
    auto a32 = Fixed64::Max();
    auto a33 = Fixed64::MaxInt();
    auto a34 = Fixed64::Zero();
    auto a35 = Fixed64::One();

    TestNaNInf(count);
    TestFromInt(count);
    TestFromFloat(count);
    TestFromFraction(count);
    TestToInt(count);
    TestToFloorInt(count);
    TestToCeilInt(count);
    TestToRoundInt(count);
    TestAdd(count);
    TestSub(count);
    TestMul(count);
    TestDiv(count);
}