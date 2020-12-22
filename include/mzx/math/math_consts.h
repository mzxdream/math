#ifndef __MZX_MATH_CONSTS_H__
#define __MZX_MATH_CONSTS_H__

#include <algorithm>

namespace mzx
{
    template <typename T, typename = void>
    class MathConsts;

    template <>
    class MathConsts<float>
    {
    public:
        static constexpr float EPSILON = 1e-06f;
        static constexpr float PI = 3.1415926535897932384626433832795f;
        static constexpr float HALF_PI = PI * 0.5f;
        static constexpr float TWO_PI = PI * 2.0f;
        static constexpr float SQRT2 = 1.4142135623730950488016887242097f;
        static constexpr float HALF_SQRT2 = SQRT2 * 0.5f;
        static constexpr float TWO_SQRT2 = SQRT2 * 2.0f;
    };

    template <>
    class MathConsts<double>
    {
    public:
        static constexpr double EPSILON = 1e-08;
        static constexpr double PI = 3.1415926535897932384626433832795;
        static constexpr double HALF_PI = PI * 0.5;
        static constexpr double TWO_PI = PI * 2.0;
        static constexpr double SQRT2 = 1.4142135623730950488016887242097;
        static constexpr double HALF_SQRT2 = SQRT2 * 0.5;
        static constexpr double TWO_SQRT2 = SQRT2 * 2.0;
    };

    template <typename T>
    class MathConsts<T, std::enable_if_t<std::is_integral_v<T>>>
    {
    public:
        using RType = T;

    public:
        static constexpr RType EPSILON = 0;
        static constexpr RType PI = 3;
        static constexpr RType HALF_PI = PI / 2;
        static constexpr RType TWO_PI = PI * 2;
        static constexpr RType SQRT2 = 1;
        static constexpr RType HALF_SQRT2 = SQRT2 / 2;
        static constexpr RType TWO_SQRT2 = SQRT2 * 2;
    };
} // namespace mzx

#endif