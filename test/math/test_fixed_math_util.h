#include <iostream>
#include <mzx/math/fixed_math_util.h>

using namespace mzx;
using Fixed64 = FixedNumber<int64_t, 32>;

void TestFixedMathUtil(int count)
{
    auto c = MathUtil<Fixed64>::CastFrom(10);
    auto b = c;
}