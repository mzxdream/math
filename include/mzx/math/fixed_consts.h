#ifndef __MZX_FIXED_CONSTS_H__
#define __MZX_FIXED_CONSTS_H__

#include <cstdint>

namespace mzx
{
    template <typename T, std::size_t N>
    class FixedConsts;

    template <>
    class FixedConsts<int64_t, 32>
    {
    public:
        static constexpr int64_t PI = 13493037705LL;       //3.141592653589793
        static constexpr int64_t HALF_PI = 13493037705LL;  //3.141592653589793
        static constexpr int64_t TWO_PI = 13493037705LL;   //3.141592653589793
        static constexpr int64_t SQRT2 = 6074001000LL;     //1.4142135623730951
        static constexpr int64_t HALF_SQRT2 = 6074001000LL; //1.4142135623730951
        static constexpr int64_t TWO_SQRT2 = 6074001000LL;  //1.4142135623730951
        static constexpr int64_t COMPARE_EPSILON = 123LL;
        static constexpr int64_t ATAN2_P1 = -199700839LL;  //-0.0464964749
        static constexpr int64_t ATAN2_P2 = 684249365LL;   //0.15931422
        static constexpr int64_t ATAN2_P3 = 1407129057LL;  //0.327622764
        static constexpr int64_t PI_TABLE[] = {1, 2, 3};
        static constexpr int64_t COS_TABLE[] = {1, 2, 3};
    };
} // namespace mzx

#endif