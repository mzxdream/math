#ifndef __MZX_FIXED_NUMBER_H__
#define __MZX_FIXED_NUMBER_H__

#include <cstdint>
#include <cmath>
#include <limits>
#include <cassert>
#include <mzx/math/fixed_consts.h>

namespace mzx
{
    template <typename T, std::size_t N, typename F>
    class FixedNumber
    {
    public:
        static_assert(N > 0 && N < sizeof(T));
        using RType = std::enable_if_t<std::is_integral_v<T> && std::is_signed_v<T>, T>;
        using RUType = std::make_unsigned_t<RType>;
        using FType = std::enable_if_t<std::is_floating_point_v<F>, F>;

    public:
        static constexpr RType R_NBITS = N;
        static constexpr RType R_BASE = 1 << R_NBITS;
        static constexpr RType R_MASK = R_BASE - 1;
        static constexpr RType R_HALF = R_BASE >> 1;
        static constexpr RType R_NAN = std::numeric_limits<RType>::min();
        static constexpr RType R_INF = std::numeric_limits<RType>::max();
        static constexpr RType R_MAX = R_INF - 1;
        static constexpr RType R_MAX_INT = R_MAX >> R_NBITS;
        static constexpr FType R_MAX_FLT = static_cast<FType>(R_MAX) / static_cast<FType>(R_BASE);

    public:
        explicit FixedNumber(RType raw_value)
            : raw_value_(raw_value)
#ifndef NDEBUG
              ,
              debug_value_(RawToFloat(raw_value))
#endif
        {
        }

    public:
        void Set(RType raw_value)
        {
            raw_value_ = raw_value;
#ifndef NDEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
        }
        RType Get() const
        {
            return raw_value_;
        }
        FType ToFloat() const
        {
            return RawToFloat(raw_value_);
        }
        RType ToInt() const
        {
            assert(raw_value_ >= -R_MAX && raw_value_ <= R_MAX);
            return raw_value_ >= 0 ? (raw_value_ >> R_NBITS) : -(-raw_value_ >> R_NBITS);
        }
        RType ToFloorInt() const
        {
            assert(raw_value_ >= -R_MAX && raw_value_ <= R_MAX);
            return raw_value_ >> FBITS;
        }
        RType ToCeilInt() const
        {
            assert(raw_value_ >= -R_MAX && raw_value_ <= R_MAX);
            if (raw_value_ >= 0)
            {
                return static_cast<RType>((static_cast<RUType>(raw_value_) + R_MASK) >> R_NBITS);
            }
            return -(-raw_value_ >> R_NBITS);
        }
        RType ToRoundInt() const
        {
            assert(raw_value_ >= -R_MAX && raw_value_ <= R_MAX);
            if (raw_value_ >= 0)
            {
                return static_cast<RType>((static_cast<RUType>(raw_value_) + R_HALF) >> R_NBITS);
            }
            return -static_cast<RType>((static_cast<RUType>(-raw_value_) + R_HALF) >> R_NBITS);
        }
        FixedNumber ToFloor() const
        {
            return FixedNumber::FromInt(ToFloorInt());
        }
        FixedNumber ToCeil() const
        {
            return FixedNumber::FromInt(ToCeilInt());
        }
        FixedNumber ToRound() const
        {
            return FixedNumber::FromInt(ToRoundInt());
        }
        bool IsNAN() const
        {
            return raw_value_ == R_NAN;
        }
        bool IsPosInfinite() const
        {
            return raw_value_ == R_INF;
        }
        bool IsNegInfinite() const
        {
            return raw_value_ == -R_INF;
        }
        bool IsInfinite() const
        {
            return raw_value_ == R_INF || raw_value_ == -R_INF;
        }
        bool IsFinite() const
        {
            return raw_value_ != R_NAN && raw_value_ != R_INF && raw_value_ != -R_INF;
        }

    public:
        FixedNumber operator+(const FixedNumber &a) const
        {
            return FixedNumber(FixedAdd(raw_value_, a.raw_value_));
        }
        FixedNumber operator-(const FixedNumber &a) const
        {
            return FixedNumber(FixedSub(raw_value_, a.raw_value_));
        }
        FixedNumber operator-() const
        {
            return IsNAN() ? *this : FixedNumber(-raw_value_);
        }
        FixedNumber operator*(const FixedNumber &a) const
        {
            return FixedNumber(FixedMul(raw_value_, a.raw_value_));
        }
        FixedNumber operator/(const FixedNumber &a) const
        {
            return FixedNumber(FixedDiv(raw_value_, a.raw_value_));
        }
        FixedNumber &operator+=(const FixedNumber &a)
        {
            raw_value_ = FixedAdd(raw_value_, a.raw_value_);
#ifndef NDEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        FixedNumber &operator-=(const FixedNumber &a)
        {
            raw_value_ = FixedSub(raw_value_, a.raw_value_);
#ifndef NDEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        FixedNumber &operator*=(const FixedNumber &a)
        {
            raw_value_ = FixedMul(raw_value_, a.raw_value_);
#ifndef NDEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        FixedNumber &operator/=(const FixedNumber &a)
        {
            raw_value_ = FixedDiv(raw_value_, a.raw_value_);
#ifndef NDEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        bool operator==(const FixedNumber &a) const
        {
            return raw_value_ != R_NAN && a.raw_value_ != R_NAN && raw_value_ == a.raw_value_;
        }
        bool operator!=(const FixedNumber &a) const
        {
            return !(*this == a);
        }
        bool operator<(const FixedNumber &a) const
        {
            return raw_value_ != R_NAN && a.raw_value_ != R_NAN && raw_value_ < a.raw_value_;
        }
        bool operator<=(const FixedNumber &a) const
        {
            return raw_value_ != R_NAN && a.raw_value_ != R_NAN && raw_value_ <= a.raw_value_;
        }
        bool operator>(const FixedNumber &a) const
        {
            return raw_value_ != R_NAN && a.raw_value_ != R_NAN && raw_value_ > a.raw_value_;
        }
        bool operator>=(const FixedNumber &a) const
        {
            return raw_value_ != R_NAN && a.raw_value_ != R_NAN && raw_value_ >= a.raw_value_;
        }

    public:
        static FixedNumber FromInt(RType value)
        {
            assert(value >= -R_MAX_INT && value <= R_MAX_INT);
            return FixedNumber(value << R_BITS);
        }
        static FixedNumber FromFloat(FType value)
        {
            switch (std::fpclassify(value))
            {
            case FP_NAN:
                return FixedNumber(R_NAN);
            case FP_INFINITE:
                return value < 0 ? FixedNumber(-R_INF) : FixedNumber(R_INF);
            default:
                assert(value >= -R_MAX_FLT && value <= -R_MAX_FLT);
                return FixedNumber(static_cast<RType>(value * R_BASE));
            }
        }
        static FixedNumber FromFraction(RType numerator, RType denominator)
        {
            assert(numerator >= -R_MAX_INT && numerator <= R_MAX_INT && denominator >= -R_MAX_INT && denominator <= R_MAX_INT && denominator != 0);
            RUType p = std::abs(numerator);
            RUType q = std::abs(denominator);
            auto r = static_cast<RType>(((p << (R_BITS + 1)) / q + 1) >> 1);
            return FixedNumber((numerator ^ denominator) < 0 ? -r : r);
        }
        static const FixedNumber &Nan()
        {
            static const FixedNumber a(R_NAN);
            return a;
        }
        static const FixedNumber &Infinity()
        {
            static const FixedNumber a(R_INF);
            return a;
        }
        static const FixedNumber &Epsilon()
        {
            static const FixedNumber a(1);
            return a;
        }
        static const FixedNumber &Min()
        {
            static const FixedNumber a(1);
            return a;
        }
        static const FixedNumber &Max()
        {
            static const FixedNumber a(R_MAX);
            return a;
        }
        static const FixedNumber &MaxInt()
        {
            static const FixedNumber a(R_MAX & ~R_MASK);
            return a;
        }
        static const FixedNumber &Zero()
        {
            static const FixedNumber a(0);
            return a;
        }
        static const FixedNumber &One()
        {
            static const FixedNumber a(R_BASE);
            return a;
        }

    private:
        static FType RawToFloat(RType raw_value)
        {
            switch (raw_value)
            {
            case R_NAN:
                return std::numeric_limits<FType>::quiet_NaN();
            case -R_INF:
                return -std::numeric_limits<FType>::infinity();
            case R_INF:
                return std::numeric_limits<FType>::infinity();
            default:
                return static_cast<FType>(raw_value) / static_cast<FType>(R_BASE);
            }
        }
        static RType FixedAdd(RType a, RType b)
        {
            if (a == R_NAN || b == R_NAN)
            {
                return R_NAN;
            }
            if (a == R_INF || a == -R_INF)
            {
                if (a == -b)
                {
                    return R_NAN;
                }
                return a;
            }
            if (b == R_INF || b == -R_INF)
            {
                return b;
            }
            assert(a + b - a == b && a + b >= -R_MAX && a + b <= R_MAX);
            return a + b;
        }
        static int64_t FixedSub(int64_t a, int64_t b)
        {
            if (a == FNAN || b == FNAN)
            {
                return FNAN;
            }
            if (a == FINF || a == -FINF)
            {
                if (a == b)
                {
                    return FNAN;
                }
                return a;
            }
            if (b == FINF || b == -FINF)
            {
                return -b;
            }
            MZX_CHECK(std::abs(static_cast<double>(a) - b) <= FMAX);
            return a - b;
        }
        static int64_t FixedMul(int64_t a, int64_t b)
        {
            if (a == FNAN || b == FNAN)
            {
                return FNAN;
            }
            uint64_t p = std::abs(a);
            uint64_t q = std::abs(b);
            if (p == FINF || q == FINF)
            {
                if (a == 0 || b == 0)
                {
                    return FNAN;
                }
                return (a ^ b) < 0 ? -FINF : FINF;
            }
            uint64_t x1 = p >> FBITS;
            uint64_t y1 = q >> FBITS;
            uint64_t x2 = p & FMASK;
            uint64_t y2 = q & FMASK;
            MZX_CHECK((double)x1 * y1 * FBASE <= FMAX);
            MZX_CHECK(((double)x2 * y2 + FHALF) / 2 <= FMAX);
            MZX_CHECK((((double)x1 * y1) * FBASE + (double)x1 * y2 + (double)x2 * y1 + (((double)x2 * y2 + FHALF) / FBASE)) <= FMAX);
            int64_t r = static_cast<int64_t>(((x1 * y1) << FBITS) + x1 * y2 + x2 * y1 + ((x2 * y2 + FHALF) >> FBITS));
            return (a ^ b) < 0 ? -r : r;
        }
        static int64_t FixedDiv(int64_t a, int64_t b)
        {
            if (a == FNAN || b == FNAN)
            {
                return FNAN;
            }
            uint64_t dividend = std::abs(a);
            uint64_t divisor = std::abs(b);
            if (dividend == FINF)
            {
                if (divisor == FINF)
                {
                    return FNAN;
                }
                return (a ^ b) < 0 ? -FINF : FINF;
            }
            if (divisor == FINF)
            {
                return 0;
            }
            if (b == 0)
            {
                return a == 0 ? FNAN : (a > 0 ? FINF : -FINF);
            }
            uint64_t result = 0;
            int64_t bits = FBITS + 1;
            while ((divisor & 0xF) == 0 && bits >= 4)
            {
                divisor >>= 4;
                bits -= 4;
            }
            while (dividend != 0 && bits >= 0)
            {
                while ((dividend & 0xF000000000000000ULL) == 0 && bits >= 4)
                {
                    dividend <<= 4;
                    bits -= 4;
                }
                while ((dividend & 0x8000000000000000ULL) == 0 && bits >= 1)
                {
                    dividend <<= 1;
                    bits--;
                }
                uint64_t t = dividend / divisor;
                dividend %= divisor;
                MZX_CHECK(((double)t * (1ULL << bits)) / 2 <= FMAX);
                MZX_CHECK(((double)result + (t << bits)) / 2 <= FMAX);
                result += (t << bits);
                dividend <<= 1;
                --bits;
            }
            return (a ^ b) < 0 ? -static_cast<int64_t>((result + 1) >> 1) : static_cast<int64_t>((result + 1) >> 1);
        }

    private:
        int64_t raw_value_;
#ifdef MZX_DEBUG
        float debug_value_;
#endif
    };
} // namespace mzx

#endif