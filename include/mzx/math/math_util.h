#ifndef __MZX_MATH_UTIL_H__
#define __MZX_MATH_UTIL_H__

#include <cmath>
#include <mzx/math/math_consts.h>

namespace mzx
{
    template <typename T, typename = void>
    class MathUtil;

    template <typename T>
    class MathUtil<T, std::enable_if_t<std::is_integral_v<T>>>
    {
        using RType = T;

    public:
        static constexpr RType PI()
        {
            return 3;
        }
        static constexpr RType HalfPI()
        {
            return 2;
        }
        static constexpr RType TwoPI()
        {
            return 6;
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

    template <typename T>
    class MathUtil<T, std::enable_if_t<std::is_floating_point_v<T>>>
    {
        using RType = T;

    public:
        static constexpr RType PI()
        {
            return 3.1415926535897932384626433832795;
        }
        static constexpr RType HalfPI()
        {
            return PI() * 0.5;
        }
        static constexpr RType TwoPI()
        {
            return PI() * 2.0;
        }
        static constexpr RType Sqrt2()
        {
            return 1.4142135623730950488016887242097;
        }
        static constexpr RType HalfSqrt2()
        {
            return Sqrt2() / 2;
        }
        static constexpr RType TwoSqrt2()
        {
            return Sqrt2() * 2.0;
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
            return std::abs(a - b) < 1e-06 ? 0 : (a < b ? -1 : 1);
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

    using MathUtilF = MathUtil<float>;
    using MathUtilD = MathUtil<double>;
    using MathUtilI = MathUtil<int>;
} // namespace mzx

#endif