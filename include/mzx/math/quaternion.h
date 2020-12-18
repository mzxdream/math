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
        void SetX(const RType &x)
        {
            x_ = x;
        }
        const RType &Y() const
        {
            return y_;
        }
        void SetY(const RType &y)
        {
            y_ = y;
        }
        const RType &Z() const
        {
            return z_;
        }
        void SetZ(const RType &z)
        {
            z_ = z;
        }
        const RType &W() const
        {
            return w_;
        }
        void SetW(const RType &w)
        {
            w_ = w;
        }
        void Set(const RType &x, const RType &y, const RType &z, const RType &w)
        {
            x_ = x;
            y_ = y;
            z_ = z;
            w_ = w;
        }
        void Set(const RType *arr)
        {
            x_ = arr[0];
            y_ = arr[1];
            z_ = arr[2];
            w_ = arr[3];
        }

    public:
        RType &operator[](int i)
        {
            return (&x_)[i];
        }
        const RType &operator[](int i) const
        {
            return &x_[i];
        }

    public:
        static const Quaternion &Identity()
        {
            static const Quaternion a(R_ZERO, R_ZERO, R_ZERO, R_ONE);
            return a;
        }

    private:
        RType x_;
        RType y_;
        RType z_;
        RType w_;
    };
} // namespace mzx

#endif