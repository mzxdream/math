#ifndef __MZX_FIXED_MATH_UTIL_H__
#define __MZX_FIXED_MATH_UTIL_H__

#include <mzx/math/fixed_consts.h>
#include <mzx/math/fixed_number.h>
#include <mzx/math/math_util.h>

namespace mzx
{
    template <typename T, std::size_t N, typename F>
    class MathUtil<FixedNumber<T, N, F>>
    {
        using RType = FixedNumber<T, N, F>;
        using RConsts = FixedConsts<T, N>;

    public:
        static const RType &PI()
        {
            static const RType a(RConsts::PI);
            return a;
        }
        static const RType &HalfPI()
        {
            static const RType a(RConsts::HALF_PI);
            return a;
        }
        static const RType &TwoPI()
        {
            static const RType a(RConsts::TWO_PI);
            return a;
        }
        static constexpr RType Sqrt2()
        {
            return 1;
        }
        static constexpr RType HalfSqrt2()
        {
            return 1;
        }
        static constexpr RType TwoSqrt2()
        {
            return 3;
        }
        static constexpr RType CastFrom(int a)
        {
            return static_cast<RType>(a);
        }
        static constexpr RType CastFrom(int numerator, int denominator)
        {
            return static_cast<RType>(numerator) / static_cast<RType>(denominator);
        }
        static constexpr int CompareApproximately(RType a, RType b)
        {
            return a == b ? 0 : (a < b ? -1 : 1);
        }
        static constexpr RType Abs(RType a)
        {
            return std::abs(a);
        }
        static constexpr RType Min(RType a, RType b)
        {
            return std::min(a, b);
        }
        static constexpr RType Max(RType a, RType b)
        {
            return std::max(a, b);
        }
        static constexpr RType Clamp(RType a, RType mina, RType maxa)
        {
            return std::min(std::max(a, mina), maxa);
        }
        static constexpr RType Lerp(RType a, RType b, RType t)
        {
            return a + (b - a) * t;
        }
        static constexpr RType Sqrt(RType a)
        {
            return sqrt(a);
        }
        static constexpr RType Cos(RType a)
        {
            return cos(a);
        }
        static constexpr RType Sin(RType a)
        {
            return sin(a);
        }
        static constexpr RType Acos(RType a)
        {
            return acos(a);
        }
        static constexpr RType Asin(RType a)
        {
            return asin(a);
        }
        static constexpr RType Atan(RType a)
        {
            return atan(a);
        }
        static constexpr RType Atan2(RType b, RType a)
        {
            return atan2(b, a);
        }
        static constexpr RType Rad2Deg(RType rad)
        {
            return rad * static_cast<RType>(180) / PI();
        }
        static constexpr RType Deg2Rad(RType deg)
        {
            return deg * PI() / static_cast<RType>(180);
        }
    };
} // namespace mzx

#endif