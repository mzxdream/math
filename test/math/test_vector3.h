#include <iostream>
#include <mzx/math/vector3.h>

using namespace mzx;

inline void TestVector3(int count)
{
    Vector3F a;
    Vector3F b(1.0f, 1.2f, 1.0f);
    Vector3F c(2.0f, 2.0f, 3.0f);

    a.SetX(1);
    a.SetY(2);
    a.SetZ(3);
    a.Set(0, 0, 0);
    a.Set(3, 3, 3);
    a.Scale(c);
    a.Normalize();
    a.Set(10, 11, 12);
    auto k1 = a.SqrMagnitude();
    auto k2 = a.Magnitude();
    a.Set(0, 0, 0);

    auto k3 = a[0];
    a[1] = 2;
    //
    auto t1 = a + b;
    auto t2 = a - b;
    auto t3 = a * 2.0f;
    auto t4 = 2.0f * a;
    auto t5 = a / 2.0f;
    a += b;
    a -= b;
    a *= 2.0f;
    a /= 2.0f;
    auto c1 = a == b;
    auto c2 = a != b;

    a = Vector3F(1, 1, 1);
    b = Vector3F(15, 15, 15);
    auto p1 = Vector3F::Lerp(a, b, 0.5f);
    auto p2 = Vector3F::LerpUnclamped(a, b, 2.0f);
    auto p3 = Vector3F::MoveTowards(a, b, 3.0f);
    auto cur_vel = Vector3F::One();
    auto stime = 1.0f;
    auto speed = 1.0f;
    auto delta_time = 1.0f;
    auto p4 = Vector3F::SmoothDamp(a, b, cur_vel, stime, speed, delta_time);
    auto p5 = Vector3F::Scale(a, b);
    auto p6 = Vector3F::Cross(a, b);
    auto p7 = Vector3F::Reflect(a, b);
    auto p8 = Vector3F::Normalize(b);
    auto p9 = Vector3F::Dot(a, b);
    auto p10 = Vector3F::Project(a, b);
    auto p11 = Vector3F::ProjectOnPlane(a, b);
    auto p12 = Vector3F::Angle(a, b);
    auto p13 = Vector3F::SignedAngle(a, b, Vector3F::Up());
    auto p14 = Vector3F::SqrDistance(a, b);
    auto p15 = Vector3F::Distance(a, b);
    auto p16 = Vector3F::ClampMagnitude(a, 10.0f);
    auto p17 = Vector3F::SqrMagnitude(a);
    auto p18 = Vector3F::Magnitude(a);
    auto p19 = Vector3F::Min(a, b);
    auto p20 = Vector3F::Max(a, b);
    auto p21 = Vector3F::Slerp(a, b, 2.0f);
    auto p22 = Vector3F::SlerpUnclamped(a, b, 2.0f);
    Vector3F::OrthoNormalize(&a, &b);
    Vector3F::OrthoNormalize(&a, &b, &c);
    auto p23 = Vector3F::RotateTowards(a, b, 2.0f, 2.0f);

    auto b1 = Vector3F::Zero();
    auto b2 = Vector3F::One();
    auto b3 = Vector3F::Forward();
    auto b4 = Vector3F::Back();
    auto b5 = Vector3F::Up();
    auto b6 = Vector3F::Down();
    auto b7 = Vector3F::Left();
    auto b8 = Vector3F::Right();
    auto b9 = b7 + b8;
}