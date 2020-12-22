#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include <mzx/math/math_util.h>

using namespace mzx;

inline void TestMathUtil(int count)
{
   auto c = MathUtil<float>::Abs(1.0f);
   auto d = MathUtil<float>::PI;
   auto  e = d;
}