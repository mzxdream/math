#include <iostream>
#include <mzx/math/quaternion.h>

using namespace mzx;

inline void TestQuaternion(int count)
{
    QuaternionF q1;
    QuaternionF q2(1, 2, 3, 4);
    q1.SetLookRotation(Vector3F(2, 0, 0));
    auto a1 = q1.EulerAngles();
    q1.SetEulerAngles(Vector3F(30, 30, 30));
    float angle = 0;
    Vector3F axis;
    q1.ToAngleAxis(angle, axis);
    q1.SetFromToRotation(Vector3F(1, 0, 0), Vector3F(-1, 0, 0));
    q1.Normalize();
    auto a2 = q1.Normalized();
    auto c1 = q1[0];
    auto c2 = q1[2];
    q1[2] = 10;
    auto a3 = q1 * q2;
    auto a4 = q1 * Vector3F(1, 1, 1);
    auto a5 = q1 == q2;
    auto a6 = q1 != q2;
    auto a7 = QuaternionF::Dot(q1, q2);
    auto a8 = QuaternionF::Angle(q1, q2);
    auto a9 = QuaternionF::Euler(1, 2, 3);
    auto a10 = QuaternionF::Euler(Vector3F(1, 2, 3));
    auto a11 = QuaternionF::RotateTowards(q1, q2, 2.0f);
    auto a12 = QuaternionF::Normalize(q1);
    Vector3F v1(2, 3, 4);
    Vector3F v2(4, 5, 6);
    auto a13 = QuaternionF::FromToRotation(v1, v2);
    auto a14 = QuaternionF::Inverse(q1);
    auto a15 = QuaternionF::Slerp(q1, q2, 0.5f);
    auto a16 = QuaternionF::Lerp(q1, q2, 0.5f);
    angle = 180.0f;
    axis = Vector3F(2, 2, 2);
    auto a17 = QuaternionF::AngleAxis(angle, axis);
    auto a18 = QuaternionF::LookRotation(axis, Vector3F::Up());
    auto a19 = QuaternionF::Identity();
}