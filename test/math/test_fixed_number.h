#include <iostream>
#include <mzx/math/fixed_number.h>

using namespace mzx;

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
    a -=b ;
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
}