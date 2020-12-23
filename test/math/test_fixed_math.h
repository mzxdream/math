#include <iostream>
#include <mzx/math/fixed_math.h>

using namespace mzx;

void TestFixedMath(int count)
{
    auto c = MathUtil<Fixed64>::CastFrom(10);
    auto b = c;
}