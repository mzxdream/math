#include <iostream>
#include <mzx/math/fixed_math_util.h>

using namespace mzx;

void TestFixedMath(int count)
{
    auto c = MathUtil<Fixed64>::CastFrom(10);
    auto b = c;
}