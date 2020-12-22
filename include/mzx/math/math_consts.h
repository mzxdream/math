#ifndef __MZX_MATH_CONSTS_H__
#define __MZX_MATH_CONSTS_H__

#include <algorithm>

namespace mzx
{
    template <typename T, typename = void>
    class MathConsts;

    template <>
    class MathConsts<float>
    {
    public:
        static constexpr float Epsilon()
        {
            return 1e-06f;
        }
        static constexpr float PI()
        {
            return 3.1415926535897932384626433832795f;
        }
        static constexpr float HalfPI()
        {
            return PI() * 0.5f;
        }
        static constexpr float TwoPI()
        {
            return PI() * 2.0f;
        }
        static constexpr float Sqrt2()
        {
            return 1.4142135623730950488016887242097f;
        }
        static constexpr float HalfSqrt2()
        {
            return Sqrt2() * 0.5f;
        }
        static constexpr float TwoSqrt2()
        {
            return Sqrt2() * 2.0f;
        }
    };

    template <>
    class MathConsts<double>
    {
    public:
        static constexpr double Epsilon()
        {
            return 1e-08;
        }
        static constexpr double PI()
        {
            return 3.1415926535897932384626433832795;
        }
        static constexpr double HalfPI()
        {
            return PI() * 0.5;
        }
        static constexpr double TwoPI()
        {
            return PI() * 2.0;
        }
        static constexpr double Sqrt2()
        {
            return 1.4142135623730950488016887242097;
        }
        static constexpr double HalfSqrt2()
        {
            return Sqrt2() * 0.5;
        }
        static constexpr double TwoSqrt2()
        {
            return Sqrt2() * 2.0;
        }
    };

    template <typename T>
    class MathConsts<T, std::enable_if_t<std::is_integral_v<T>>>
    {
    public:
        using RType = T;

    public:
        static constexpr RType Epsilon()
        {
            return 0;
        }
        static constexpr RType PI()
        {
            return 3;
        }
        static constexpr RType HalfPI()
        {
            return 1;
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
    };

    using MathConstsF = MathConsts<float>;
    using MathConstsD = MathConsts<double>;
    using MathConstsI = MathConsts<int>;
} // namespace mzx

#endif