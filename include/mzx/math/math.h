#ifndef __MZX_MATH_H__
#define __MZX_MATH_H__

#include <cstdint>

namespace mzx
{
    template <typename T>
    inline T Abs(T a)
    {
        return std::abs(a);
    }
} // namespace mzx

#endif