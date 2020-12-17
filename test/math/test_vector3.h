#include <iostream>
#include <mzx/math/vector3.h>

using namespace mzx;

inline void TestVector3(int count)
{
    Vector3f a;
    Vector3f b(1.0f, 1.2f, 1.0f);
    Vector3f c(2.0f, 2.0f, 3.0f);

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

    a = Vector3f(1, 1, 1);
    b = Vector3f(15, 15, 15);
    auto p1 = Vector3f::Lerp(a, b, 0.5f);
    auto p2 = Vector3f::LerpUnclamped(a, b, 2.0f);
    auto p3 = Vector3f::MoveTowards(a, b, 3.0f);
    auto cur_vel = Vector3f::One();
    auto stime = 1.0f;
    auto speed = 1.0f;
    auto delta_time = 1.0f;
    auto p4 = Vector3f::SmoothDamp(a, b, cur_vel, stime, speed, delta_time);
    auto p5 = Vector3f::Scale(a, b);
    auto p6 = Vector3f::Cross(a, b);
    auto p7 = Vector3f::Reflect(a, b);
    auto p8 = Vector3f::Normalize(b);
    auto p9 = Vector3f::Dot(a, b);
    auto p10 = Vector3f::Project(a, b);
    auto p11 = Vector3f::ProjectOnPlane(a, b);
    auto p12 = Vector3f::Angle(a, b);
    auto p13 = Vector3f::SignedAngle(a, b, Vector3f::Up());
    auto p14 = Vector3f::SqrDistance(a, b);
    auto p15 = Vector3f::Distance(a, b);
    auto p16 = Vector3f::ClampMagnitude(a, 10.0f);
    auto p17 = Vector3f::SqrMagnitude(a);
    auto p18 = Vector3f::Magnitude(a);
    auto p19 = Vector3f::Min(a, b);
    auto p20 = Vector3f::Max(a, b);

    auto b1 = Vector3f::Zero();
    auto b2 = Vector3f::One();
    auto b3 = Vector3f::Forward();
    auto b4 = Vector3f::Back();
    auto b5 = Vector3f::Up();
    auto b6 = Vector3f::Down();
    auto b7 = Vector3f::Left();
    auto b8 = Vector3f::Right();
    auto b9 = b7 + b8;
}