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
        using RCType = typename RType::RType;
        using RCUType = typename RType::RUType;
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
            RCUType bit = static_cast<RCUType>(1) << (sizeof(RCType) * CHAR_BIT - 2);
            while (bit > num)
            {
                bit >>= 2;
            }
            while (bit != 0)
            {
                if (num >= res + bit)
                {
                    num -= (res + bit);
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
            static_assert(RType::R_NBITS >= 2);
            bit = static_cast<RCUType>(1) << (RType::R_NBITS - 2);
            while (bit != 0)
            {
                if (num >= res + bit)
                {
                    num -= (res + bit);
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
        static RType CosDeg(const RType &a)
        {
            return !a.IsFinite() ? RType::Nan() : RType(CosDegRaw(a.Get()));
        }
        static RType SinDeg(const RType &a)
        {
            static_assert(RType::R_BASE * 90 / 90 == RType::R_BASE);
            assert(a.Get() + RType::R_BASE * 90 - RType::R_BASE * 90 == a.Get());
            return !a.IsFinite() ? RType::Nan() : RType(-CosDegRaw(a.Get() + RType::R_BASE * 90));
        }
        static RType Cos(const RType &a)
        {
            return CosDeg(Rad2Deg(a));
        }
        static RType Sin(const RType &a)
        {
            return SinDeg(Rad2Deg(a));
        }
        static RType Acos(const RType &a)
        {
            if (Abs(a) == RType::One())
            {
                return Atan2(RType::Zero(), a);
            }
            return Atan2(Sqrt((RType::One() + a) * (RType::One() - a)), a);
        }
        static RType Asin(const RType &a)
        {
            if (Abs(a) == RType::One())
            {
                return Atan2(a, RType::Zero());
            }
            return Atan2(a, Sqrt((RType::One() + a) * (RType::One() - a)));
        }
        static RType Atan(const RType &a)
        {
            return Atan2(a, RType::One());
        }
        static RType Atan2(const RType &b, const RType &a)
        {
            if (!a.IsFinite() || !b.IsFinite())
            {
                return RType::Nan();
            }
            static const RType ATAN2_P1(RConsts::ATAN2_P1);
            static const RType ATAN2_P2(RConsts::ATAN2_P2);
            static const RType ATAN2_P3(RConsts::ATAN2_P3);
            RType x = Abs(a);
            RType y = Abs(b);
            RType t = x > y ? y / x : x / y;
            RType s = t * t;
            RType r = ((ATAN2_P1 * s + ATAN2_P2) * s - ATAN2_P3) * s * t + t;
            if (y > x)
            {
                r = HalfPI() - r;
            }
            if (a < RType::Zero())
            {
                r = PI() - r;
            }
            if (b < RType::Zero())
            {
                r = -r;
            }
            return r;
        }
        static RType Rad2Deg(const RType &rad)
        {
            if (!rad.IsFinite())
            {
                return RType::Nan();
            }
            auto raw_value = rad.Get();
            if (raw_value < 0)
            {
                for (auto i = RConsts::PI_TABLE_LEN - 1; i >= 0; --i)
                {
                    if (-raw_value >= RConsts::PI_TABLE[i])
                    {
                        raw_value += RConsts::PI_TABLE[i];
                    }
                }
            }
            else
            {
                for (auto i = RConsts::PI_TABLE_LEN - 1; i > 0; --i)
                {
                    if (raw_value > RConsts::PI_TABLE[i])
                    {
                        raw_value -= RConsts::PI_TABLE[i];
                    }
                }
            }
            static_assert(RConsts::R_BASE * 360 / 360 == RConsts::R_BASE);
            assert(raw_value * RConsts::R_BASE * 360 / (RConsts::R_BASE * 360) == raw_value);
            return RType(raw_value * RConsts::R_BASE * 360 / RConsts::TWO_PI);
        }
        static RType Deg2Rad(const RType &deg)
        {
            if (!deg.IsFinite())
            {
                return RType::Nan();
            }
            static_assert(RType::R_BASE * 360 / 360 == RType::R_BASE);
            auto raw_value = deg.Get() % (RType::R_BASE * 360);
            if (raw_value < 0)
            {
                raw_value += (RType::R_BASE * 360);
            }
            assert(raw_value * RConsts::TWO_PI / RConsts::TWO_PI == raw_value);
            return RType(raw_value * RConsts::TWO_PI / (RType::R_BASE * 360));
        }

    private:
        static RCType CosDegLookupTable(RCType deg) //[0 - 90]
        {
            static constexpr auto NUM90 = RType::R_BASE * 90;
            assert(deg >= 0 && deg <= NUM90);
            deg *= RConsts::COS_TABLE_LEN;
            auto a = deg / NUM90;
            auto b = deg - a * NUM90;
            if (b == 0)
            {
                return RConsts::COS_TABLE[a];
            }
            return (RConsts::COS_TABLE[a] * (NUM90 - b) + RConsts::COS_TABLE[a + 1] * b) / NUM90;
        }
        static RCType CosDegRaw(RCType raw_value)
        {
            if (raw_value < 0)
            {
                raw_value = -raw_value;
            }
            static_assert(RType::R_BASE * 360 / 360 == RType::R_BASE);
            if (raw_value >= RType::R_BASE * 360)
            {
                raw_value %= (RType::R_BASE * 360);
            }
            if (raw_value <= RType::R_BASE * 90)
            {
                return CosDegLookupTable(raw_value);
            }
            else if (raw_value <= RType::R_BASE * 180)
            {
                return -CosDegLookupTable(RType::R_BASE * 180 - raw_value);
            }
            else if (raw_value <= RType::R_BASE * 270)
            {
                return -CosDegLookupTable(raw_value - RType::R_BASE * 180);
            }
            return CosDegLookupTable(RType::R_BASE * 360 - raw_value);
        }
    };
} // namespace mzx

#endif