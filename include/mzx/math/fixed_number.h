#ifndef __MZX_FIXED_NUMBER_H__
#define __MZX_FIXED_NUMBER_H__

#include <cstdint>
#include <cmath>
#include <limits>
#include <cassert>
#include <mzx/math/fixed_consts.h>

namespace mzx
{
    template <typename T, std::size_t N>
    class FixedNumber
    {
        static_assert(std::is_same_v<T, std::int16_t> || std::is_same_v<T, std::int32_t> || std::is_same_v<T, std::int64_t>);
        static_assert(N > 0 && N + 1 < sizeof(T));

    public:
        using RType = T;

    public:
        static constexpr RType R_BITS = N;
        static constexpr RType R_BASE = (static_cast<RType>(1) << R_BITS);
        static constexpr RType R_MASK = R_BASE - static_cast<RType>(1);
        static constexpr RType R_HALF = R_BASE >> static_cast<RType>(1);
        static constexpr RType R_NAN = std::numeric_limits<RType>::min();
        static constexpr RType R_INF = std::numeric_limits<RType>::max();
        static constexpr RType R_MAX = R_INF - 1;
        static constexpr RType R_MAX_INT = (R_MAX >> R_BITS);
        static constexpr RType R_MAX_FLT = static_cast<float>(R_MAX) / static_cast<float>(R_BASE);

    public:
        explicit FixedNumber(int64_t raw_value)
            : raw_value_(raw_value)
#ifdef MZX_DEBUG
              ,
              debug_value_(RawToFloat(raw_value))
#endif
        {
        }

    public:
        void Set(int64_t raw_value)
        {
            raw_value_ = raw_value;
#ifdef MZX_DEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
        }
        int64_t Get() const
        {
            return raw_value_;
        }
        float ToFloat() const
        {
            return RawToFloat(raw_value_);
        }
        int64_t ToInt() const
        {
            MZX_CHECK(std::abs(raw_value_) <= FMAX);
            if (raw_value_ >= 0)
            {
                return raw_value_ >> FBITS;
            }
            return -(-raw_value_ >> FBITS);
        }
        int64_t ToFloorInt() const
        {
            MZX_CHECK(std::abs(raw_value_) <= FMAX);
            return raw_value_ >> FBITS;
        }
        int64_t ToCeilInt() const
        {
            MZX_CHECK(std::abs(raw_value_) <= FMAX);
            if (raw_value_ >= 0)
            {
                return static_cast<int64_t>((static_cast<uint64_t>(raw_value_) + FMASK) >> FBITS);
            }
            return -(-raw_value_ >> FBITS);
        }
        int64_t ToRoundInt() const
        {
            MZX_CHECK(std::abs(raw_value_) <= FMAX);
            if (raw_value_ >= 0)
            {
                return static_cast<int64_t>((static_cast<uint64_t>(raw_value_) + FHALF) >> FBITS);
            }
            return -static_cast<int64_t>((static_cast<uint64_t>(-raw_value_) + FHALF) >> FBITS);
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
            return raw_value_ == FNAN;
        }
        bool IsPosInf() const
        {
            return raw_value_ == FINF;
        }
        bool IsNegInf() const
        {
            return raw_value_ == -FINF;
        }
        bool IsInf() const
        {
            return raw_value_ == FINF || raw_value_ == -FINF;
        }
        bool IsFinite() const
        {
            return raw_value_ != FNAN && raw_value_ != FINF && raw_value_ != -FINF;
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
#ifdef MZX_DEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        FixedNumber &operator-=(const FixedNumber &a)
        {
            raw_value_ = FixedSub(raw_value_, a.raw_value_);
#ifdef MZX_DEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        FixedNumber &operator*=(const FixedNumber &a)
        {
            raw_value_ = FixedMul(raw_value_, a.raw_value_);
#ifdef MZX_DEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        FixedNumber &operator/=(const FixedNumber &a)
        {
            raw_value_ = FixedDiv(raw_value_, a.raw_value_);
#ifdef MZX_DEBUG
            debug_value_ = RawToFloat(raw_value_);
#endif
            return *this;
        }
        bool operator==(const FixedNumber &a) const
        {
            return raw_value_ != FNAN && a.raw_value_ != FNAN && raw_value_ == a.raw_value_;
        }
        bool operator!=(const FixedNumber &a) const
        {
            return raw_value_ == FNAN || a.raw_value_ == FNAN || raw_value_ != a.raw_value_;
        }
        bool operator<(const FixedNumber &a) const
        {
            return raw_value_ != FNAN && a.raw_value_ != FNAN && raw_value_ < a.raw_value_;
        }
        bool operator<=(const FixedNumber &a) const
        {
            return raw_value_ != FNAN && a.raw_value_ != FNAN && raw_value_ <= a.raw_value_;
        }
        bool operator>(const FixedNumber &a) const
        {
            return raw_value_ != FNAN && a.raw_value_ != FNAN && raw_value_ > a.raw_value_;
        }
        bool operator>=(const FixedNumber &a) const
        {
            return raw_value_ != FNAN && a.raw_value_ != FNAN && raw_value_ >= a.raw_value_;
        }

    public:
        static FixedNumber FromInt(int64_t value)
        {
            MZX_CHECK(std::abs(value) <= FMAX_INT);
            return FixedNumber(value << FBITS);
        }
        static FixedNumber FromFloat(float value)
        {
            switch (std::fpclassify(value))
            {
            case FP_NAN:
                return FixedNumber(FNAN);
            case FP_INFINITE:
                return value < 0 ? FixedNumber(-FINF) : FixedNumber(FINF);
            default:
                MZX_CHECK(std::abs(value) <= FMAX_INT);
                return FixedNumber(static_cast<int64_t>(value * FBASE));
            }
        }
        static FixedNumber FromFraction(int64_t numerator, int64_t denominator)
        {
            MZX_CHECK(std::abs(numerator) <= FMAX_INT && denominator != 0);
            uint64_t p = std::abs(numerator);
            uint64_t q = std::abs(denominator);
            int64_t r = static_cast<int64_t>(((p << (FBITS + 1)) / q + 1) >> 1);
            return FixedNumber((numerator ^ denominator) < 0 ? -r : r);
        }
        static const FixedNumber &NaN()
        {
            static const FixedNumber a(FNAN);
            return a;
        }
        static const FixedNumber &Infinity()
        {
            static const FixedNumber a(FINF);
            return a;
        }
        static const FixedNumber &Epsilon()
        {
            static const FixedNumber epsilon(1);
            return epsilon;
        }
        static const FixedNumber &Min()
        {
            static const FixedNumber a(1);
            return a;
        }
        static const FixedNumber &Max()
        {
            static const FixedNumber a(FMAX);
            return a;
        }
        static const FixedNumber &MaxInt()
        {
            static const FixedNumber a(FMAX_INT);
            return a;
        }
        static const FixedNumber &Zero()
        {
            static const FixedNumber a(0);
            return a;
        }
        static const FixedNumber &One()
        {
            static const FixedNumber a(FBASE);
            return a;
        }

    private:
        static float RawToFloat(int64_t raw_value)
        {
            switch (raw_value)
            {
            case -FINF:
                return -std::numeric_limits<float>::infinity();
            case FINF:
                return std::numeric_limits<float>::infinity();
            case FNAN:
                return std::numeric_limits<float>::quiet_NaN();
            default:
                return static_cast<float>(raw_value) / FBASE;
            }
        }
        static int64_t FixedAdd(int64_t a, int64_t b)
        {
            if (a == FNAN || b == FNAN)
            {
                return FNAN;
            }
            if (a == FINF || a == -FINF)
            {
                if (a == -b)
                {
                    return FNAN;
                }
                return a;
            }
            if (b == FINF || b == -FINF)
            {
                return b;
            }
            MZX_CHECK(std::abs(static_cast<double>(a) + b) <= FMAX);
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