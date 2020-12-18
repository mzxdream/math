#ifndef __MZX_QUATERNION_H__
#define __MZX_QUATERNION_H__

#include <cmath>
#include <mzx/math/vector3.h>

namespace mzx
{
    template <typename T>
    class QuaternionMathUtil
    {
    public:
        using RType = T;

    public:
        static RType Epsilon();
        static RType Acos(const RType &a);
        static RType Rad2Deg(const RType &rad);
        static RType Deg2Rad(const RType &deg);
    };

    template <typename T>
    class Quaternion
    {
    public:
        using RType = T;
        using MathUtil = QuaternionMathUtil<RType>;

    private:
        static const RType R_EPSILON; //1e-6f
        static const RType R_SQR_EPSILON;
        static const RType R_FLIP; //1e-4f
        static const RType R_ZERO; //0
        static const RType R_ONE;  //1
        static const RType R_TWO;  //2
        static const RType R_360;  //360

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
        void SetLookRotation(const Vector3<RType> &view, const Vector3<RType> &up = Vector3<RType>::Up())
        {
            *this = LookRotation(view, up);
        }
        Vector3<RType> EulerAngles() const
        {
            auto vec3 = ToEulerRad(*this);
            vec3.Set(
                MakePosAngle(MathUtil::Rad2Deg(vec3.X())),
                MakePosAngle(MathUtil::Rad2Deg(vec3.Y())),
                MakePosAngle(MathUtil::Rad2Deg(vec3.Z())));
            return vec3;
        }
        void SetEulerAngles(const Vector3<RType> &a)
        {
            *this = FromEulerRad(Vector3<RType>(
                MathUtil::Deg2Rad(a.X()),
                MathUtil::Deg2Rad(a.Y()),
                MathUtil::Deg2Rad(a.Z())));
        }
        void ToAngleAxis(RType *angle, Vector3<RType> *axis)
        {
            ToAngleAxisRad(*this, angle, axis);
            if (angle != nullptr)
            {
                *angle = MathUtil::Rad2Deg(*angle);
            }
        }
        void SetFromToRotation(const Vector3 &from, const Vector3 &to)
        {
            *this = FromToRotation(from, to);
        }
        void Normalize()
        {
            *this = Normalize(*this);
        }
        Quaternion Normalized() const
        {
            return Normalize(*this);
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
        Quaternion operator*(const Quaternion &a) const
        {
            return Quaternion(
                w_ * a.x_ + x_ * a.w_ + y_ * a.z_ - z_ * a.y_,
                w_ * a.y_ + y_ * a.w_ + z_ * a.x_ - x_ * a.z_,
                w_ * a.z_ + z_ * a.w_ + x_ * a.y_ - y_ * a.x_,
                w_ * a.w_ - x_ * a.x_ - y_ * a.y_ - z_ * a.z_);
        }
        Quaternion operator*(const Vector3<RType> &point) const
        {
            auto x = x_ * R_TWO;
            auto y = y_ * R_TWO;
            auto z = z_ * R_TWO;
            auto xx = x_ * x;
            auto yy = y_ * y;
            auto zz = z_ * z;
            auto xy = x_ * y;
            auto xz = x_ * z;
            auto yz = y_ * z;
            auto wx = w_ * x;
            auto wy = w_ * y;
            auto wz = w_ * z;

            return Vector3<RType>(
                (R_ONE - (yy + zz)) * point.X() + (xy - wz) * point.Y() + (xz + wy) * point.Z(),
                (xy + wz) * point.X() + (R_ONE - (xx + zz)) * point.Y() + (yz - wx) * point.Z(),
                (xz - wy) * point.X() + (yz + wx) * point.Y() + (R_ONE - (xx + yy)) * point.Z());
        }
        bool operator==(const Quaternion &a) const
        {
            return IsEqualUsingDot(Dot(*this, a));
        }
        bool operator!=(const Quaternion &a) const
        {
            return !(*this == a);
        }

    public:
        static RType Dot(const Quaternion &a, const Quaternion &b)
        {
            return a.x_ * b.x_ + a.y_ * b.y_ + a.z_ * b.z_ + a.w_ * b.w_;
        }
        static RType Angle(const Quaternion &a, const Quaternion &b)
        {
            auto dot = Dot(a, b);
            if (IsEqualUsingDot(dot))
            {
                return R_ZERO;
            }
            return MathUtil::Rad2Deg(MathUtil::Acos(RMin(RAbs(dot), R_ONE)) * R_TWO);
        }
        static Quaternion Euler(const RType &x, const RType &y, const RType &z)
        {
            return FromEulerRad(Vector3<RType>(
                MathUtil::Deg2Rad(x),
                MathUtil::Deg2Rad(y),
                MathUtil::Deg2Rad(z)));
        }
        static Quaternion Euler(const Vector3<RType> &euler)
        {
            return Euler(euler.X(), euler.Y(), euler.Z());
        }
        static Quaternion RotateTowards(const Quaternion &from, const Quaternion &to, const RType &max_degrees)
        {
            auto angle = Angle(from, to);
            if (angle <= R_EPSILON)
            {
                return to;
            }
            return SlerpUnclamped(from, to, RMin(R_ONE, max_degrees / angle));
        }
        static Quaternion Normalize(const Quaternion &a)
        {
            auto mag = MathUtil::Sqrt(Dot(a, a));
            if (mag <= R_EPSILON)
            {
                return Identity();
            }
            return Quaternion(a.x_ / mag, a.y_ / mag, a.z_ / mag, a.w_ / mag);
        }

        static Quaternion FromToRotation(const Vector3<RType> &from, const Vector3<RType> &to)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion Inverse(const Quaternion &a)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion Slerp(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion SlerpUnclamped(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion Lerp(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion LerpUnclamped(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion AngleAxis(const RType &angle, const Vector3<RType> &axis)
        {
            //TODO
            return Quaternion();
        }
        static Quaternion LookRotation(const Vector3<RType> &view, const Vector3<RType> &up)
        {
            //TODO
            return Quaternion();
        }
        static const Quaternion &Identity()
        {
            static const Quaternion a(R_ZERO, R_ZERO, R_ZERO, R_ONE);
            return a;
        }

    private:
        static bool IsEqualUsingDot(const RType &dot)
        {
            return dot > R_ONE - R_EPSILON;
        }
        static Quaternion FromEulerRad(const Vector3<RType> &a)
        {
            //TODO
            return Quaternion();
        }
        static Vector3<RType> ToEulerRad(const Quaternion &a);
        {
            //TODO
            return Vector3<RType>();
        }
        static void ToAngleAxisRad(const Quaternion &a, RType *angle, Vector3<RType> *axis)
        {
            //TODO
        }
        static RType MakePosAngle(const RType &a)
        {
            static const RType NEG_FLIP = MathUtil::Rad2Deg(-R_FLIP);
            static const RType POS_FLIP = R_360 - R_FLIP;
            if (a < NEG_FLIP)
            {
                return a + R_360;
            }
            if (a > POS_FLIP)
            {
                return a - R_360;
            }
            return a;
        }
        static RType RAbs(const RType &a)
        {
            return a >= R_ZERO ? a : -a;
        }
        static RType RMin(const RType &a, const RType &b)
        {
            return a < b ? a : b;
        }
        static RType RMax(const RType &a, const RType &b)
        {
            return a > b ? a : b;
        }

    private:
        RType x_;
        RType y_;
        RType z_;
        RType w_;
    };
} // namespace mzx

#endif