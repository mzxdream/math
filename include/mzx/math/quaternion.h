#ifndef __MZX_QUATERNION_H__
#define __MZX_QUATERNION_H__

#include <cmath>

namespace mzx
{
    template <typename T>
    class QuaternionMathUtil
    {
    public:
        using RType = T;

    public:
        static RType Epsilon();
    };

    template <typename T>
    class Quaternion
    {
    public:
        using RType = T;
        using MathUtil = QuaternionMathUtil<RType>;

    private:
        static const RType R_EPSILON;
        static const RType R_SQR_EPSILON;
        static const RType R_ZERO;
        static const RType R_ONE;

    public:
        Quaternion()
            : x_(R_ZERO), y_(R_ZERO), z_(R_ZERO), w_(R_ZERO)
        {
        }
        explicit Quaternion(const RType &x, const RType &y, const RType &z, const RType &w)
            : x_(x), y_(y), z_(z), w_(w)
        {
        }

    public:
        const RType &X() const
        {
            return x_;
        }

    public:
    private:
        RType x_;
        RType y_;
        RType z_;
        RType w_;
    };
} // namespace mzx

#endif