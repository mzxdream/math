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
        using RCType = RType::RType;
        using RCUType = RType::RUType;
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
        static const RType &Sqrt2()
        {
            static const RType a(RConsts::SQRT2);
            return a;
        }
        static const RType &HalfSqrt2()
        {
            static const RType a(RConsts::HALF_SQRT2);
            return a;
        }
        static const RType &TwoSqrt2()
        {
            static const RType a(RConsts::TWO_SQRT2);
            return a;
        }
        static RType CastFrom(int a)
        {
            return RType::FromInt(a);
        }
        static RType CastFrom(int numerator, int denominator)
        {
            return RType::FromFraction(numerator, denominator);
        }
        static int CompareApproximately(const RType &a, const RType &b)
        {
            static const RType t(RConsts::COMPARE_EPSILON);
            return Abs(a - b) < t ? 0 : (a < b ? -1 : 1);
        }
        static RType Abs(const RType &a)
        {
            if (a.IsNAN())
            {
                return a;
            }
            auto mask = a.Get() >> (sizeof(a.Get()) * CHAR_BIT - 1);
            return RType((a.Get() + mask) ^ mask);
        }
        static RType Min(const RType &a, const RType &b)
        {
            return a < b ? a : b;
        }
        static RType Max(const RType &a, const RType &b)
        {
            return a < b ? b : a;
        }
        static RType Clamp(const RType &a, const RType &mina, const RType &maxa)
        {
            return a < mina ? mina : (a > maxa ? maxa : a);
        }
        static RType Lerp(const RType &a, const RType &b, const RType &t)
        {
            return a + (b - a) * t;
        }
        static RType Sqrt(const RType &a)
        {
            if (a.IsNAN() || a.IsPosInf())
            {
                return a;
            }
            if (a < RType::Zero())
            {
                return RType::Nan();
            }
            if (a == RType::Zero())
            {
                return RType::Zero();
            }
            RCUType num = a.Get();
            RCUType res = 0;
            RCUType bit = static_cast<RCUType>(1) << (sizeof(RCUType) * CHAR_BIT - 2);
            while (bit > num)
            {
                bit >>= 2;
            }
            while (bit != 0)
            {
                if (num >= res + bit)
                {
                    num -= res + bit;
                    res = (res >> 1) + bit;
                }
                else
                {
                    res >>= 1;
                }
                bit >>= 2;
            }
            if (num > RType::R_MASK)
            {
                num -= res;
                num = (num << RType::R_NBITS) - RType::R_HALF;
                res = (res << RType::R_NBITS) + RType::R_HALF;
            }
            else
            {
                num <<= RType::R_NBITS;
                res <<= RType::R_NBITS;
            }
            bit = static_cast<RCUType>(1) << (RType::R_NBITS - 2);
            while (bit != 0)
            {
                if (num >= res + bit)
                {
                    num -= res + bit;
                    res = (res >> 1) + bit;
                }
                else
                {
                    res >>= 1;
                }
                bit >>= 2;
            }
            if (num > res)
            {
                ++res;
            }
            return RType(static_cast<RCType>(res));
        }
        static RType Pow(const RType &a, int n)
        {
            if (n == 0)
            {
                return RType::One();
            }
            if (n == 1 || a == RType::Zero() || a == RType::One())
            {
                return a;
            }
            if (n < 0)
            {
                return RType::One() / Pow(a, -n);
            }
            return (n & 1) == 0 ? Pow(a * a, n >> 1) : Pow(a * a, n >> 1) * a;
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