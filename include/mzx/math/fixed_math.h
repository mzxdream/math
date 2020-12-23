#ifndef __MZX_FIXED_MATH_H__
#define __MZX_FIXED_MATH_H__

#include <mzx/math/fixed_number.h>
#include <mzx/math/math_util.h>

namespace mzx
{
    template <typename T, std::size_t N, typename F>
    class MathUtil<FixedNumber<T, N, F>>
    {
        using RType = FixedNumber<T, N, F>;

    public:
        static const RType CastFrom(int a)
        {
            return RType::FromInt(a);
        }
    };
} // namespace mzx

#endif