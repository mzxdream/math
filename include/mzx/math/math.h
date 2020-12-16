#ifndef __MZX_MATH_H__
#define __MZX_MATH_H__

#include <cstdint>
#include <algorithm>

namespace mzx
{
    template <typename T>
    inline T Abs(T a)
    {
        return std::abs(a);
    }
} // namespace mzx

#endif