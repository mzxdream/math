#ifndef __MZX_VECTOR3_H__
#define __MZX_VECTOR3_H__

#include <cmath>

namespace mzx
{
    template <typename T>
    class Vector3MathUtil
    {
    public:
        using RType = T;

    public:
        static RType PI();
        static RType Epsilon();
        static RType Zero();
        static RType One();
        static RType Two();
        static RType SmoothTimeMin();
        static RType Dot48();
        static RType Dot235();
        static RType Sqrt(const RType &a);
        static RType Rad2Deg(const RType &a);
        static RType Deg2Rad(const RType &a);
        static RType Acos(const RType &a);
    };

    template <typename T>
    class Vector3
    {
    public:
        using RType = T;
        using MathUtil = Vector3MathUtil<RType>;

    private:
        static const RType R_PI;
        static const RType R_EPSILON;
        static const RType R_SQR_EPSILON;
        static const RType R_ZERO;
        static const RType R_ONE;
        static const RType R_TWO;
        static const RType R_SMOOTH_TIME_MIN;
        static const RType R_DOT48;
        static const RType R_DOT235;

    public:
        Vector3()
            : x_(R_ZERO), y_(R_ZERO), z_(R_ZERO)
        {
        }
        explicit Vector3(const RType &x, const RType &y, const RType &z)
            : x_(x), y_(y), z_(z)
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
        void Set(const RType &x, const RType &y, const RType &z)
        {
            x_ = x;
            y_ = y;
            z_ = z;
        }
        void Set(const RType *arr)
        {
            x_ = arr[0];
            y_ = arr[1];
            z_ = arr[2];
        }
        void Scale(const Vector3 &s)
        {
            x_ *= s.x_;
            y_ *= s.y_;
            z_ *= s.z_;
        }
        void Normalize()
        {
            *this = Normalize(*this);
        }
        Vector3 Normalized() const
        {
            return Normalize(*this);
        }
        RType SqrMagnitude() const
        {
            return x_ * x_ + y_ * y_ + z_ * z_;
        }
        RType Magnitude() const
        {
            return MathUtil::Sqrt(x_ * x_ + y_ * y_ + z_ * z_);
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
        Vector3 operator+(const Vector3 &a) const
        {
            return Vector3(x_ + a.x_, y_ + a.y_, z_ + a.z_);
        }
        Vector3 operator-(const Vector3 &a) const
        {
            return Vector3(x_ - a.x_, y_ - a.y_, z_ - a.z_);
        }
        Vector3 operator-() const
        {
            return Vector3(-x_, -y_, -z_);
        }
        Vector3 operator*(const RType &s) const
        {
            return Vector3(x_ * s, y_ * s, z_ * s);
        }
        friend Vector3 operator*(const RType &s, const Vector3 &a)
        {
            return Vector3(a.x_ * s, a.y_ * s, a.z_ * s);
        }
        Vector3 operator/(const RType &d) const
        {
            return Vector3(x_ / d, y_ / d, z_ / d);
        }
        Vector3 &operator+=(const Vector3 &a)
        {
            x_ += a.x_;
            y_ += a.y_;
            z_ += a.z_;
            return *this;
        }
        Vector3 &operator-=(const Vector3 &a)
        {
            x_ -= a.x_;
            y_ -= a.y_;
            z_ -= a.z_;
            return *this;
        }
        Vector3 &operator*=(const RType &s)
        {
            x_ *= s;
            y_ *= s;
            z_ *= s;
            return *this;
        }
        Vector3 &operator/=(const RType &d)
        {
            x_ /= d;
            y_ /= d;
            z_ /= d;
            return *this;
        }
        bool operator==(const Vector3 &a) const
        {
            return SqrDistance(*this, a) <= R_SQR_EPSILON;
        }
        bool operator!=(const Vector3 &a) const
        {
            return !(*this == a);
        }

    public:
        static Vector3 Lerp(const Vector3 &a, const Vector3 &b, const RType &t)
        {
            return LerpUnclamped(a, b, RClamp(t, R_ZERO, R_ONE));
        }
        static Vector3 LerpUnclamped(const Vector3 &a, const Vector3 &b, const RType &t)
        {
            return Vector3(a.x_ + (b.x_ - a.x_) * t, a.y_ + (b.y_ - a.y_) * t, a.z_ + (b.z_ - a.z_) * t);
        }
        static Vector3 MoveTowards(const Vector3 &current, const Vector3 &target, const RType &max_distance)
        {
            auto x = target.x_ - current.x_;
            auto y = target.y_ - current.y_;
            auto z = target.z_ - current.z_;
            auto sqr_dist = x * x + y * y + z * z;
            if (sqr_dist <= R_SQR_EPSILON || (max_distance >= R_ZERO && sqr_dist <= max_distance * max_distance))
            {
                return target;
            }
            auto dist = MathUtil::Sqrt(sqr_dist);
            return Vector3(current.x_ + x * max_distance / dist, current.y_ + y * max_distance / dist, current.z_ + z * max_distance / dist);
        }
        static Vector3 SmoothDamp(const Vector3 &current, const Vector3 &target, Vector3 &current_velocity, const RType &smooth_time, const RType &max_speed, const RType &delta_time)
        {
            auto stime = RMax(smooth_time, R_SMOOTH_TIME_MIN);
            auto omega = R_TWO / stime;
            auto x = omega * delta_time;
            auto exp = R_ONE / (R_ONE + x + R_DOT48 * x * x + R_DOT235 * x * x * x);

            auto change = current - target;
            auto max_change = max_speed * stime;
            auto sqr_mag = change.SqrMagnitude();
            if (sqr_mag > max_change * max_change)
            {
                change = change * max_change / MathUtil::Sqrt(sqr_mag);
            }

            auto temp = (current_velocity + change * omega) * delta_time;
            current_velocity = (current_velocity - temp * omega) * exp;
            auto output = (current - change) + (change + temp) * exp;

            auto orig_minus_current = target - current;
            auto out_minus_orig = output - target;
            if (Dot(orig_minus_current, out_minus_orig) > R_ZERO)
            {
                output = target;
                current_velocity = (output - target) / delta_time;
            }
            return output;
        }
        static Vector3 Scale(const Vector3 &a, const Vector3 &b)
        {
            return Vector3(a.x_ * b.x_, a.y_ * b.y_, a.z_ * b.z_);
        }
        static Vector3 Cross(const Vector3 &lhs, const Vector3 &rhs)
        {
            return Vector3(lhs.y_ * rhs.z_ - lhs.z_ * rhs.y_, lhs.z_ * rhs.x_ - lhs.x_ * rhs.z_, lhs.x_ * rhs.y_ - lhs.y_ * rhs.x_);
        }
        static Vector3 Reflect(const Vector3 &in_direction, const Vector3 &in_normal)
        {
            auto factor = -R_TWO * Dot(in_normal, in_direction);
            return Vector3(factor * in_normal.x_ + in_direction.x_, factor * in_normal.y_ + in_direction.y_, factor * in_normal.z_ + in_direction.z_);
        }
        static Vector3 Normalize(const Vector3 &value)
        {
            auto mag = Magnitude(value);
            if (mag > R_EPSILON)
            {
                return value / mag;
            }
            return Vector3(R_ZERO, R_ZERO, R_ZERO);
        }
        static RType Dot(const Vector3 &lhs, const Vector3 &rhs)
        {
            return lhs.x_ * rhs.x_ + lhs.y_ * rhs.y_ + lhs.z_ * rhs.z_;
        }
        static Vector3 Project(const Vector3 &vec, const Vector3 &on_normal)
        {
            auto sqr_mag = on_normal.SqrMagnitude();
            if (sqr_mag > R_EPSILON)
            {
                return on_normal * Dot(vec, on_normal) / sqr_mag;
            }
            return Vector3(R_ZERO, R_ZERO, R_ZERO);
        }
        static Vector3 ProjectOnPlane(const Vector3 &vec, const Vector3 &plane_normal)
        {
            auto sqr_mag = plane_normal.SqrMagnitude();
            if (sqr_mag > R_EPSILON)
            {
                return vec - plane_normal * Dot(vec, plane_normal) / sqr_mag;
            }
            return vec;
        }
        static RType Angle(const Vector3 &from, const Vector3 &to)
        {
            auto denominator = MathUtil::Sqrt(from.SqrMagnitude() * to.SqrMagnitude());
            if (denominator > R_SQR_EPSILON)
            {
                auto dot = RClamp(Dot(from, to) / denominator, -R_ONE, R_ONE);
                return MathUtil::Rad2Deg(MathUtil::Acos(dot));
            }
            return R_ZERO;
        }
        static RType SignedAngle(const Vector3 &from, const Vector3 &to, const Vector3 &axis)
        {
            auto cross_x = from.y_ * to.z_ - from.z_ * to.y_;
            auto cross_y = from.z_ * to.x_ - from.x_ * to.z_;
            auto cross_z = from.x_ * to.y_ - from.y_ * to.x_;
            if (axis.x_ * cross_x + axis.y_ * cross_y + axis.z_ * cross_z >= R_ZERO)
            {
                return Angle(from, to);
            }
            return -Angle(from, to);
        }
        static RType SqrDistance(const Vector3 &a, const Vector3 &b)
        {
            return (a - b).SqrMagnitude();
        }
        static RType Distance(const Vector3 &a, const Vector3 &b)
        {
            return (a - b).Magnitude();
        }
        static Vector3 ClampMagnitude(const Vector3 &vec, const RType &max_magnitude)
        {
            auto sqr_mag = vec.SqrMagnitude();
            if (sqr_mag > max_magnitude * max_magnitude)
            {
                auto mag = MathUtil::Sqrt(sqr_mag);
                return Vector3(vec.x_ * max_magnitude / mag, vec.y_ * max_magnitude / mag, vec.z_ * max_magnitude / mag);
            }
            return vec;
        }
        static RType SqrMagnitude(const Vector3 &vec)
        {
            return vec.x_ * vec.x_ + vec.y_ * vec.y_ + vec.z_ * vec.z_;
        }
        static RType Magnitude(const Vector3 &vec)
        {
            return MathUtil::Sqrt(vec.x_ * vec.x_ + vec.y_ * vec.y_ + vec.z_ * vec.z_);
        }
        static Vector3 Min(const Vector3 &lhs, const Vector3 &rhs)
        {
            return Vector3(RMin(lhs.x_, rhs.x_), RMin(lhs.y_, rhs.y_), RMin(lhs.z_, rhs.z_));
        }
        static Vector3 Max(const Vector3 &lhs, const Vector3 &rhs)
        {
            return Vector3(RMax(lhs.x_, rhs.x_), RMax(lhs.y_, rhs.y_), RMax(lhs.z_, rhs.z_));
        }
        static Vector3 Slerp(const Vector3 &lhs, const Vector3 &rhs, const RType &t)
        {
            return SlerpUnclamped(lhs, rhs, RClamp(t, R_ZERO, R_ONE));
        }
        static Vector3 SlerpUnclamped(const Vector3 &lhs, const Vector3 &rhs, const RType &t)
        {
            auto lhs_mag = Magnitude(lhs);
            auto rhs_mag = Magnitude(rhs);
            if (lhs_mag < R_EPSILON || rhs_mag < R_EPSILON)
            {
                return Lerp(lhs, rhs, t);
            }
            auto lerped_magnitude = Lerp(lhs_mag, rhs_mag, t);
            auto dot = Dot(lhs, rhs) / (lhs_mag * rhs_mag);
            if (dot > R_ONE - R_EPSILON)
            {
                return Lerp(lhs, rhs, t);
            }
            else if (dot < -R_ONE + R_EPSILON)
            {
                auto lhs_norm = lhs / lhs_mag;
                auto axis = OrthoNormalVectorFast(lhs_norm);
                RType m[3][3];
                m.SetAxisAngle(axis, R_PI * t);
                auto slerped = m.MultiplyPoint3(lhs_norm);
                slerped *= lerped_magnitude;
                return slerped;
            }
            else
            {
                auto axis = Cross(lhs, rhs);
                auto lhs_norm = lhs / lhs_mag;
                axis = Normalize(axis);
                auto angle = MathUtil::Acos(dot) * t;
                RType m[3][3];
                m.SetAxisAngle(axis, angle);
                auto slerped = m.MultiplyPoint3(lhs_norm);
                slerped *= lerped_magnitude;
                return slerped;
            }
        }
        static void OrthoNormalize(Vector3 &normal, Vector3 &tangent)
        {
        }
        static void OrthoNormalize(Vector3 &normal, Vector3 &tangent, Vector3 &binormal)
        {
        }
        static Vector3 RotateTowards(const Vector3 &current, const Vector3 &target, const RType &max_radians, const RType &max_magnitude)
        {
        }
        static const Vector3 &Zero()
        {
            static const Vector3 a(R_ZERO, R_ZERO, R_ZERO);
            return a;
        }
        static const Vector3 &One()
        {
            static const Vector3 a(R_ONE, R_ONE, R_ONE);
            return a;
        }
        static const Vector3 &Forward()
        {
            static const Vector3 a(R_ZERO, R_ZERO, R_ONE);
            return a;
        }
        static const Vector3 &Back()
        {
            static const Vector3 a(R_ZERO, R_ZERO, -R_ONE);
            return a;
        }
        static const Vector3 &Up()
        {
            static const Vector3 a(R_ZERO, R_ONE, R_ZERO);
            return a;
        }
        static const Vector3 &Down()
        {
            static const Vector3 a(R_ZERO, -R_ONE, R_ZERO);
            return a;
        }
        static const Vector3 &Left()
        {
            static const Vector3 a(-R_ONE, R_ZERO, R_ZERO);
            return a;
        }
        static const Vector3 &Right()
        {
            static const Vector3 a(R_ONE, R_ZERO, R_ZERO);
            return a;
        }

    private:
        static RType RMin(const RType &a, const RType &b)
        {
            return a < b ? a : b;
        }
        static RType RMax(const RType &a, const RType &b)
        {
            return a > b ? a : b;
        }
        static RType RClamp(const RType &t, const RType &mint, const RType &maxt)
        {
            return t < mint ? mint : (t > maxt ? maxt : t);
        }

    private:
        RType x_;
        RType y_;
        RType z_;
    };

    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_PI = Vector3<T>::MathUtil::PI();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_EPSILON = Vector3<T>::MathUtil::Epsilon();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_SQR_EPSILON = Vector3<T>::MathUtil::Epsilon() * Vector3<T>::MathUtil::Epsilon();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_ZERO = Vector3<T>::MathUtil::Zero();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_ONE = Vector3<T>::MathUtil::One();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_TWO = Vector3<T>::MathUtil::Two();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_SMOOTH_TIME_MIN = Vector3<T>::MathUtil::SmoothTimeMin();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_DOT48 = Vector3<T>::MathUtil::Dot48();
    template <typename T>
    const typename Vector3<T>::RType Vector3<T>::R_DOT235 = Vector3<T>::MathUtil::Dot235();

    template <>
    class Vector3MathUtil<float>
    {
    public:
        static float PI()
        {
            return 3.14159265358979323846f;
        }
        static float Epsilon()
        {
            return 1e-6f;
        }
        static float Zero()
        {
            return 0;
        }
        static float One()
        {
            return 1;
        }
        static float Two()
        {
            return 2;
        }
        static float SmoothTimeMin()
        {
            return 1e-4f;
        }
        static float Dot48()
        {
            return 0.48f;
        }
        static float Dot235()
        {
            return 0.235f;
        }
        static float Sqrt(const float &a)
        {
            return static_cast<float>(sqrt(a));
        }
        static float Rad2Deg(const float &a)
        {
            return a * 180.0f / PI();
        }
        static float Deg2Rad(const float &a)
        {
            return a * PI() / 180.0f;
        }
        static float Acos(const float &a)
        {
            return static_cast<float>(acos(a));
        }
    };

    using Vector3F = Vector3<float>;
} // namespace mzx

#endif