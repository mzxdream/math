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
        static constexpr auto PI = 13493037705LL;       //3.141592653589793
        static constexpr auto SQRT2 = 6074001000LL;     //1.4142135623730951
        static constexpr auto RAD2DEG = 246083499208LL; //57.29577951308232
        static constexpr auto DEG2RAD = 74961321LL;     //0.017453292519943295
        static constexpr auto ATAN2_P1 = -199700839LL;  //-0.0464964749
        static constexpr auto ATAN2_P2 = 684249365LL;   //0.15931422
        static constexpr auto ATAN2_P3 = 1407129057LL;  //0.327622764
    };
} // namespace mzx

#endif