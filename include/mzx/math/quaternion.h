#ifndef __MZX_QUATERNION_H__
#define __MZX_QUATERNION_H__

#include <cassert>
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
        static RType CompareApproximately(const RType &a, const RType &b);
        static RType Acos(const RType &a);
        static RType Rad2Deg(const RType &rad);
        static RType Deg2Rad(const RType &deg);
    };

    template <>
    class QuaternionMathUtil<float>
    {
    private:
        static constexpr float EPSILON = 0.000001f;
        static constexpr float MATH_PI = 3.1415926535897931f;
        static constexpr float MATH_E = 2.7182818284590451f;

    public:
        static float Epsilon()
        {
            return EPSILON;
        }
        static float PI()
        {
            return MATH_PI;
        }
        static float ConvertTo(int a)
        {
            return static_cast<float>(a);
        }
        static float ConvertTo(int numerator, int denominator)
        {
            return static_cast<float>(numerator) / denominator;
        }
        static float CompareApproximately(const float &a, const float &b)
        {
            return abs(a - b) <= EPSILON;
        }
        static float Abs(const float &a)
        {
            return abs(a);
        }
        static float Clamp(const float &a, const float &mina, const float &maxa)
        {
            return std::min(std::max(a, mina), maxa);
        }
        static float Min(const float &a, const float &b)
        {
            return std::min(a, b);
        }
        static float Max(const float &a, const float &b)
        {
            return std::max(a, b);
        }
        static float Sqrt(const float &a)
        {
            return sqrt(a);
        }
        static float Cos(const float &a)
        {
            return cos(a);
        }
        static float Sin(const float &a)
        {
            return sin(a);
        }
        static float Acos(const float &a)
        {
            return acos(a);
        }
        static float Asin(const float &a)
        {
            return asin(a);
        }
        static float Atan(const float &a)
        {
            return atan(a);
        }
        static float Atan2(const float &b, const float &a)
        {
            return atan2(b, a);
        }
        static float Rad2Deg(const float &rad)
        {
            return rad * 180.0f / MATH_PI;
        }
        static float Deg2Rad(const float &deg)
        {
            return deg * MATH_PI / 180.0f;
        }
    };

    template <typename T>
    class Quaternion
    {
    public:
        using RType = T;
        using MathUtil = QuaternionMathUtil<RType>;

    private:
        static const RType R_EPSILON;
        static const RType R_ZERO;               //0
        static const RType R_ONE;                //1
        static const RType R_TWO;                //2
        static const RType R_360;                //360
        static const RType R_DOT95;              //0.95
        static const RType R_FLIP;               //1e-4f
        static const RType R_SINGULARITY_CUTOFF; //0.499999f

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
            assert(arr != nullptr);
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
                MakeAnglePositive(MathUtil::Rad2Deg(vec3.X())),
                MakeAnglePositive(MathUtil::Rad2Deg(vec3.Y())),
                MakeAnglePositive(MathUtil::Rad2Deg(vec3.Z())));
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
            ToAxisAngleRad(*this, axis, angle);
            if (angle != nullptr)
            {
                *angle = MathUtil::Rad2Deg(*angle);
            }
        }
        void SetFromToRotation(const Vector3<RType> &from, const Vector3<RType> &to)
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
            return Dot(*this, a) > R_ONE - R_EPSILON;
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
            if (dot > R_ONE - R_EPSILON)
            {
                return R_ZERO;
            }
            return MathUtil::Rad2Deg(MathUtil::Acos(MathUtil::Min(MathUtil::Abs(dot), R_ONE)) * R_TWO);
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
            if (angle < R_EPSILON)
            {
                return to;
            }
            return SlerpUnclamped(from, to, MathUtil::Min(R_ONE, max_degrees / angle));
        }
        static Quaternion Normalize(const Quaternion &a)
        {
            auto mag = MathUtil::Sqrt(Dot(a, a));
            if (mag < R_EPSILON)
            {
                return Identity();
            }
            return Quaternion(a.x_ / mag, a.y_ / mag, a.z_ / mag, a.w_ / mag);
        }
        static Quaternion FromToRotation(const Vector3<RType> &from, const Vector3<RType> &to)
        {
            auto from_mag = from.Magnitude();
            auto to_mag = to.Magnitude();
            if (from_mag < R_EPSILON || to_mag < R_EPSILON)
            {
                return Identity();
            }
            RType m[3][3];
            SetMatrix3x3FromToRotation(m, from / from_mag, to / to_mag);
            Quaternion q;
            Matrix3x3ToQuaternion(m, q);
            return q;
        }
        static Quaternion Inverse(const Quaternion &a)
        {
            return Quaternion(-a.x_, -a.y_, -a.z_, a.w_);
        }
        static Quaternion Slerp(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            return SlerpUnclamped(a, b, MathUtil::Clamp(t, R_ZERO, R_ONE));
        }
        static Quaternion SlerpUnclamped(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            Quaternion tmp = b;
            auto dot = Dot(a, b);
            if (dot < R_ZERO)
            {
                dot = -dot;
                tmp.Set(-tmp.x_, -tmp.y_, -tmp.z_, -tmp.w_);
            }
            if (dot < R_DOT95)
            {
                auto angle = MathUtil::Acos(dot);
                auto sina = MathUtil::Sin(angle);
                auto sinat = MathUtil::Sin(angle * t);
                auto sinaomt = MathUtil::Sin(angle * (R_ONE - t));
                tmp.Set(
                    (a.x_ * sinaomt + tmp.x_ * sinat) / sina,
                    (a.y_ * sinaomt + tmp.y_ * sinat) / sina,
                    (a.z_ * sinaomt + tmp.z_ * sinat) / sina,
                    (a.w_ * sinaomt + tmp.w_ * sinat) / sina);
                return tmp;
            }
            return LerpUnclamped(a, tmp, t);
        }
        static Quaternion Lerp(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            return LerpUnclamped(a, b, MathUtil::Clamp(t, R_ZERO, R_ONE));
        }
        static Quaternion LerpUnclamped(const Quaternion &a, const Quaternion &b, const RType &t)
        {
            Quaternion tmp;
            if (Dot(a, b) < R_ZERO)
            {
                tmp.Set(
                    a.x_ + t * (-b.x_ - a.x_),
                    a.y_ + t * (-b.y_ - a.y_),
                    a.z_ + t * (-b.z_ - a.z_),
                    a.w_ + t * (-b.w_ - a.w_)));
            }
            else
            {
                tmp.Set(
                    a.x_ + t * (b.x - a.x_),
                    a.y_ + t * (b.y - a.y_),
                    a.z_ + t * (b.z - a.z_),
                    a.w_ + t * (b.w - a.w_));
            }
            return Normalize(tmp);
        }
        static Quaternion AngleAxis(const RType &angle, const Vector3<RType> &axis)
        {
            auto angle_rad = MathUtil::Deg2Rad(angle);
            auto mag = axis.Magnitude();
            if (mag > R_EPSILON)
            {
                auto a = angle / R_TWO;
                auto cosa = MathUtil::Cos(a);
                auto sina = MathUtil::Sin(a);
                return Quaternion(
                    sina * axis.X() / mag,
                    sina * axis.Y() / mag,
                    sina * axis.Z() / mag,
                    cosa);
            }
            return Identity();
        }
        static Quaternion LookRotation(const Vector3<RType> &view, const Vector3<RType> &up)
        {
            auto q = Indentity();
            if (!LookRotationToQuaternion(view, up, &q))
            {
                auto mag = view.Magnitude();
                if (mag > R_EPSILON)
                {
                    RType m[3][3];
                    SetMatrix3x3FromToRotation(m, Vector3<RType>::Forward(), view / mag);
                    Matrix3x3ToQuaternion(m, q);
                }
            }
            return q;
        }
        static const Quaternion &Identity()
        {
            static const Quaternion a(R_ZERO, R_ZERO, R_ZERO, R_ONE);
            return a;
        }

    private:
        static RType MakeAnglePositive(const RType &a)
        {
            static const RType NEG_FLIP = MathUtil::Rad2Deg(-R_FLIP);
            static const RType POS_FLIP = R_360 + NEG_FLIP;
            return a < NEG_FLIP ? (a + R_360) : (a > POS_FLIP ? (a - R_360) : a);
        }
        static Quaternion FromEulerRad(const Vector3<RType> &a)
        {
            auto x = a.x_ / R_TWO;
            auto y = a.y_ / R_TWO;
            auto z = a.z_ / R_TWO;

            auto cosx = MathUtil::Cos(x);
            auto sinx = MathUtil::Sin(x);
            auto cosy = MathUtil::Cos(y);
            auto siny = MathUtil::Sin(y);
            auto cosz = MathUtil::Cos(z);
            auto sinz = MathUtil::Sin(z);

            Quaternion qx(sinx, R_ZERO, R_ZERO, cosx);
            Quaternion qy(R_ZERO, siny, R_ZERO, cosy);
            Quaternion qz(R_ZERO, R_ZERO, sinz, cosz);

            auto ret = (qy * qx) * qz; //zxy
            assert(MathUtil::CompareApproximately(ret.SqrMagnitude(), R_ONE));
            return ret;
        }
        static Vector3<RType> ToEulerRad(const Quaternion &a)
        {
            auto q = a.Normalized();

            enum VIndexs = {X1 = 0, X2, Y1, Y2, Z1, Z2, SINGULARITY_TEST};
            enum QIndexs = {XX = 0, XY, XZ, XW, YY, YZ, YW, ZZ, ZW, WW};
            RType v[7] = {R_ZERO};
            RType d[10] = {q.x_ * q.x_, q.x_ * q.y_, q.x_ * q.z_, q.x_ * q.w_, q.y_ * q.y_, q.y_ * q.z_, q.y_ * q.w_, q.z_ * q.z_, q.z_ * q.w_, q.w_ * q.w_};
            // zxy
            v[SINGULARITY_TEST] = d[YZ] - d[XW];
            v[Z1] = R_TWO * (d[XY] + d[ZW]);
            v[Z2] = d[YY] - d[ZZ] - d[XX] + d[WW];
            v[X1] = -R_ONE;
            v[X2] = R_TWO * v[SINGULARITY_TEST];
            if (MathUtil::Abs(v[SINGULARITY_TEST]) < R_SINGULARITY_CUTOFF)
            {
                v[Y1] = TWO * (d[XZ] + d[YW]);
                v[Y2] = d[ZZ] - d[XX] - d[YY] + d[WW];
                return Vector3<RType>(
                    v[X1] * MathUtil::Asin(MathUtil::Clamp(v[X2], -R_ONE, R_ONE)),
                    MathUtil::Atan2(v[Y1], v[Y2]),
                    MathUtil::Atan2(v[Z1], v[Z2]));
            }
            auto a = d[XY] + d[ZW];
            auto b = -d[YZ] + d[XW];
            auto c = d[XY] - d[ZW];
            auto e = d[YZ] + d[XW];

            v[Y1] = a * e + b * c;
            v[Y2] = b * e - a * c;
            return Vector3<RType>(
                v[X1] * MathUtil::Asin(MathUtil::Clamp(v[X2], -R_ONE, R_ONE)),
                MathUtil::Atan2(v[Y1], v[Y2]),
                R_ZERO);
        }
        static void ToAxisAngleRad(const Quaternion &a, Vector3<RType> *axis, RType *angle)
        {
            assert(axis != nullptr && angle != nullptr);
            auto q = a.Normalized();

            *angle = R_TWO * MathUtil::Acos(q.w_);
            if (MathUtil::CompareApproximately(*angle, R_ZERO))
            {
                axis->Set(R_ONE, R_ZERO, R_ZERO);
                return;
            }
            auto t = MathUtil::Sqrt(R_ONE - q.w_ * q.w_);
            axis->Set(q.x_ / t, q.y_ / t, q.z_ / t);
        }
        static bool LookRotationToQuaternion(const Vector3<RType> &view, const Vector3<RType> &up, Quaternion *res)
        {
            assert(res != nullptr);
            RType m[3][3];
            if (!LookRotationToMatrix3x3(view, up, m))
            {
                return false;
            }
            Matrix3x3ToQuaternion(m, *res);
            return true;
        }
        static void SetMatrix3x3Indentity(RType matrix[3][3])
        {
            matrix[0][0] = R_ONE;
            matrix[0][1] = R_ZERO;
            matrix[0][2] = R_ZERO;

            matrix[1][0] = R_ZERO;
            matrix[1][1] = R_ONE;
            matrix[1][2] = R_ZERO;

            matrix[2][0] = R_ZERO;
            matrix[2][1] = R_ZERO;
            matrix[2][2] = R_ONE;
        }
        static void SetMatrix3x3Basis(RType matrix[3][3], const Vector3<RType> &inx, const Vector3<RType> &iny, const Vector3<RType> &inz)
        {
            matrix[0][0] = inx[0];
            matrix[0][1] = iny[0];
            matrix[0][2] = inz[0];

            matrix[1][0] = inx[1];
            matrix[1][1] = iny[1];
            matrix[1][2] = inz[1];

            matrix[2][0] = inx[2];
            matrix[2][1] = iny[2];
            matrix[2][2] = inz[2];
        }
        static void SetMatrix3x3FromToRotation(RType matrix[3][3], const Vector3<RType> &from, const Vector3<RType> &to)
        {
            auto v = Vector3<RType>::Cross(from, to);
            auto e = Vector3<RType>::Dot(from, to);
            if (e > R_ONE - R_EPSILON)
            {
                SetMatrix3x3Indentity(m);
            }
            else if (e < -R_ONE + R_EPSILON)
            {
                Vector3<RType> left(R_ZERO, from[2], -from[1]);
                auto dot = Vector3<RType>::Dot(left, left);
                if (dot < R_EPSILON)
                {
                    left[0] = -from[2];
                    left[1] = R_ZERO;
                    left[2] = from[0];
                }
                left /= MathUtil::Sqrt(dot);
                auto up = Vector3<RType>::Cross(left, from);

                auto fxx = -from[0] * from[0];
                auto fyy = -from[1] * from[1];
                auto fzz = -from[2] * from[2];
                auto fxy = -from[0] * from[1];
                auto fxz = -from[0] * from[2];
                auto fyz = -from[1] * from[2];

                auto uxx = up[0] * up[0];
                auto uyy = up[1] * up[1];
                auto uzz = up[2] * up[2];
                auto uxy = up[0] * up[1];
                auto uxz = up[0] * up[2];
                auto uyz = up[1] * up[2];

                auto lxx = -left[0] * left[0];
                auto lyy = -left[1] * left[1];
                auto lzz = -left[2] * left[2];
                auto lxy = -left[0] * left[1];
                auto lxz = -left[0] * left[2];
                auto lyz = -left[1] * left[2];

                matrix[0][0] = fxx + uxx + lxx;
                matrix[0][1] = fxy + uxy + lxy;
                matrix[0][2] = fxz + uxz + lxz;

                matrix[1][0] = matrix[0][1];
                matrix[1][1] = fyy + uyy + lyy;
                matrix[1][2] = fyz + uyz + lyz;

                matrix[2][0] = matrix[0][2];
                matrix[2][1] = matrix[1][2];
                matrix[2][2] = fzz + uzz + lzz;
            }
            else
            {
                auto h = (R_ONE - e) / Vector3<RType>::Dot(v, v);

                auto hvx = h * v[0];
                auto hvz = h * v[2];
                auto hvxy = hvx * v[1];
                auto hvxz = hvx * v[2];
                auto hvyz = hvz * v[1];

                matrix[0][0] = e + hvx * v[0];
                matrix[0][1] = hvxy - v[2];
                matrix[0][2] = hvxz + v[1];

                matrix[1][0] = hvxy + v[2];
                matrix[1][1] = e + h * v[1] * v[1];
                matrix[1][2] = hvyz - v[0];

                matrix[2][0] = hvxz - v[1];
                matrix[2][1] = hvyz + v[0];
                matrix[2][2] = e + hvz * v[2];
            }
        }
        static void Matrix3x3ToQuaternion(const RType matrix[3][3], Quaternion &q)
        {
            auto t = matrix[0][0] + matrix[1][1] + matrix[2][2];
            if (t > R_ZERO)
            {
                auto r = MathUtil::Sqrt(t + R_ONE);
                q.w = r / R_TWO;
                r *= R_TWO;
                q.x = (matrix[2][1] - matrix[1][2]) / r;
                q.y = (matrix[0][2] - matrix[2][0]) / r;
                q.z = (matrix[1][0] - matrix[0][1]) / r;
            }
            else
            {
                static constexpr int inext[3] = {1, 2, 0};
                int i = 0;
                if (matrix[1][1] > matrix[0][0])
                {
                    i = 1;
                }
                if (matrix[2][2] > matrix[i][i])
                {
                    i = 2;
                }
                auto j = inext[i];
                auto k = inext[j];

                auto r = MathUtil::Sqrt(matrix[i][i] - matrix[j][j] - matrix[k][k] + R_ONE);
                assert(r >= R_EPSILON);
                auto *apk_quat[3] = {&q.x_, &q.y_, &q.z_};
                *apk_quat[i] = r / TWO;
                r *= R_TWO;
                q.w = (matrix[k][j] - matrix[j][k]) / r;
                *apk_quat[j] = (matrix[j][i] + matrix[i][j]) / r;
                *apk_quat[k] = (matrix[k][i] + matrix[i][k]) / r;
            }
            q = Normalize(q);
        }
        static bool LookRotationToMatrix3x3(const Vector3<RType> &view, const Vector3<RType> &up, RType matrix[3][3])
        {
            auto z = view;

            auto mag = z.Magnitude();
            if (mag < R_EPSILON)
            {
                SetMatrix3x3Identity(matrix);
                return false;
            }
            z /= mag;

            auto x = Vector3<RType>::Cross(up, z);
            mag = x.Magnitude();
            if (mag < R_EPSILON)
            {
                SetMatrix3x3Identity(matrix);
                return false;
            }
            x /= mag;

            auto y = Vector3<RType>::Cross(z, x);
            if (!MathUtil::CompareApproximately(y.SqrMagnitude(), R_ONE))
            {
                return false;
            }
            SetMatrix3x3Basis(matrix, x, y, z);
            return true;
        }

    private:
        RType x_;
        RType y_;
        RType z_;
        RType w_;
    };
} // namespace mzx

#endif