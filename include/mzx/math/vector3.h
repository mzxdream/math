#ifndef __MZX_VECTOR3_H__
#define __MZX_VECTOR3_H__

namespace mzx
{
    template <typename T>
    class Vector3Util
    {
    public:
        static T Epsilon();
        static T EpsilonNormalSqrt();
        static T FromInt(int value);
        static T FromFraction(int numerator, int denominator);
        static T Clamp(T t, T mint, T maxt);
        static T Sqrt(T t);
    };

    template <typename T>
    class Vector3
    {
    public:
        using RType = T;

    private:
        static const RType R_EPSILON_NORMAL_SQRT = 1e-15f;
        static const RType R_EPSILON = 1e-5f;
        static const RType R_ZERO = 0;
        static const RType R_ONE = 1;
        static const RType R_TWO = 2;
        static const RType R_SMOOTH_TIME_MIN = 1e-4f;
        static const RType R_DOT48 = 0.48f;
        static const RType R_DOT235 = 0.235f;

    public:
        Vector3() = default;
        explicit Vector3(RType x, RType y, RType z)
            : x_(x), y_(y), z_(z)
        {
        }

    public:
        T X() const
        {
            return x_;
        }
        void SetX(RType x)
        {
            x_ = x;
        }
        T Y() const
        {
            return y_;
        }
        void SetY(RType y)
        {
            y_ = y;
        }
        T Z() const
        {
            return z_;
        }
        void SetZ(RType z)
        {
            z_ = z;
        }
        void Set(RType x, RType y, RType z)
        {
            x_ = x;
            y_ = y;
            z_ = z;
        }
        void Scale(const Vector3 &s)
        {
            x_ *= s.x_;
            y_ *= s.y_;
            z_ *= s.z_;
        }
        bool Equals(const Vector3 &other) const
        {
            return x_ == other.x_ && y_ == other.y_ && z_ == other.z_;
        }
        void Normalize()
        {
            auto mag = Magnitude(*this);
            if (mag > R_EPSILON)
            {
                x_ /= mag;
                y_ /= mag;
                z_ /= mag;
            }
            else
            {
                x_ = R_ZERO;
                y_ = R_ZERO;
                z_ = R_ZERO;
            }
        }
        Vector3 Normalized() const
        {
            return Normalize(*this);
        }

    public:
        RType &operator[](int i)
        {
            return &x_[i];
        }
        const RType &operator[](int i) const
        {
            return &x_[i];
        }

    public:
        static Vector3 Lerp(const Vector3 &a, const Vector3 &b, RType t)
        {
            t = Vector3Util<RType>::Clamp(t, R_ZERO, R_ONE);
            return Vector3(a.x_ + (b.x_ - a.x_) * t, a.y_ + (b.y_ - a.y_) * t, a.z_ + (b.z_ - a.z_) * t);
        }
        static Vector3 LerpUnclamped(const Vector3 &a, const Vector3 &b, RType t)
        {
            return Vector3(a.x_ + (b.x_ - a.x_) * t, a.y_ + (b.y_ - a.y_) * t, a.z_ + (b.z_ - a.z_) * t);
        }
        static Vector3 MoveTowards(const Vector3 &current, const Vector3 &target, RType max_distance)
        {
            auto x = target.x_ - current.x_;
            auto y = target.y_ - current.y_;
            auto z = target.z_ - current.z_;
            auto sqr_dist = x * x + y * y + z * z;
            if (sqr_dist == R_ZERO || (max_distance >= R_ZERO && sqr_dist <= max_distance * max_distance))
            {
                return target;
            }
            auto dist = Vector3Util<RType>::Sqrt(sqr_dist);
            return Vector3(current.x + x * max_distance / dist, current.y + y * max_distance / dist, current.z + z * max_distance / dist);
        }
        static Vector3 SmoothDamp(Vector3 current, Vector3 target, Vector3 &current_velocity, RType smooth_time, RType max_speed, RType delta_time)
        {
            if (smooth_time < R_SMOOTH_TIME_MIN)
            {
                smooth_time = R_SMOOTH_TIME_MIN;
            }

            auto omega = R_TWO / smooth_time;
            auto x = omega * delta_time;
            auto exp = R_ONE / (R_ONE + x + R_DOT48 * x * x + R_DOT235 * x * x * x);

            auto change = current - target;
            auto max_change = max_speed * smooth_time;
            auto sqr_mag = change.SqrMagnitude();
            if (sqr_mag > max_change * max_change)
            {
                change = change * max_change / Vector3Util<RType>::Sqrt(sqr_mag);
            }

            auto original_to = target;
            target = current - change;

            auto temp = (current_velocity + change * omega) * delta_time;
            current_velocity = (current_velocity - temp * omega) * exp;
            auto output = target + (change + temp) * exp;

            auto orig_minus_current = original_to - current;
            auto out_minus_orig = output - original_to;
            if (Dot(orig_minus_current, out_minus_orig) > R_ZERO)
            {
                output = original_to;
                current_velocity = (output - original_to) / delta_time;
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
            return Vector3(factor * in_normal.x + in_direction.x, factor * in_normal.y + in_direction.y, factor * in_normal.z + in_direction.z);
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

    private:
        RType x_;
        RType y_;
        RType z_;
    };
} // namespace mzx

#endif