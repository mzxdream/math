#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import math

inlContent = """#ifndef __MZX_FIXED_CONSTS_H__
#define __MZX_FIXED_CONSTS_H__

#include <cstdint>

namespace mzx
{{
    template <typename T, std::size_t N>
    class FixedConsts;

    template <>
    class FixedConsts<{rtype}, {nbits}>
    {{
    public:
        static constexpr {rtype} PI = {pi}LL; //{piRaw}
        static constexpr {rtype} HALF_PI = {halfPI}LL; //{halfPIRaw}
        static constexpr {rtype} TWO_PI = {twoPI}LL; //{twoPIRaw}
        static constexpr {rtype} SQRT2 = {sqrt2}LL; //{sqrt2Raw}
        static constexpr {rtype} HALF_SQRT2 = {halfSqrt2}LL; //{halfSqrt2Raw}
        static constexpr {rtype} TWO_SQRT2 = {twoSqrt2}LL; //{twoSqrt2Raw}
        static constexpr {rtype} COMPARE_EPSILON = {compareEpsilon}LL; //{compareEpsilonRaw}
        static constexpr {rtype} ATAN2_P1 = {atan2p1}LL; //{atan2p1Raw}
        static constexpr {rtype} ATAN2_P2 = {atan2p2}LL; //{atan2p2Raw}
        static constexpr {rtype} ATAN2_P3 = {atan2p3}LL; //{atan2p3Raw}
        //2PI,4PI,8PI,......
        static constexpr {rtype} PI_TABLE[] = {{{piTable}}};
        //[0-90]
        static constexpr {rtype} COS_TABLE[] = {{{cosTable}}};
    }};
}} // namespace mzx

#endif"""


def generateConsts(inlFilePath, rtype, nbits):
    base = (1 << nbits)

    piContent = ""
    for i in range(2, 64, 2):
        j = i * base * math.pi
        if j > 0x7FFFFFFFFFFFFFFF:
            break
        piContent += ", {0}LL".format(round(j))
    piContent = piContent[2:]

    cosContent = ""
    cosTableCount = 3600
    for i in range(cosTableCount + 1):
        cosContent += ", {0}LL".format(
            round(math.cos(i * (math.pi / 2) / cosTableCount) * base))
    cosContent = cosContent[2:]

    compareEpsilonRaw = 0.000001

    datas = {
        "rtype": rtype,
        "nbits": nbits,
        "pi": round(math.pi * base),
        "piRaw": math.pi,
        "halfPI": round(math.pi * base / 2.0),
        "halfPIRaw": math.pi / 2.0,
        "twoPI": round(math.pi * base * 2.0),
        "twoPIRaw": math.pi * 2.0,
        "sqrt2": round(math.sqrt(2) * base),
        "sqrt2Raw": math.sqrt(2),
        "halfSqrt2": round(math.sqrt(2) * base / 2.0),
        "halfSqrt2Raw": math.sqrt(2) / 2.0,
        "twoSqrt2": round(math.sqrt(2) * base * 2.0),
        "twoSqrt2Raw": math.sqrt(2) * 2.0,
        "compareEpsilon": round(compareEpsilonRaw * base),
        "compareEpsilonRaw": compareEpsilonRaw,
        "atan2p1": round(-0.0464964749 * base),
        "atan2p1Raw": -0.0464964749,
        "atan2p2": round(0.15931422 * base),
        "atan2p2Raw": 0.15931422,
        "atan2p3": round(0.327622764 * base),
        "atan2p3Raw": 0.327622764,
        "piTable": piContent,
        "cosTable": cosContent,
    }
    content = inlContent.format(**datas)
    with open(inlFilePath, "wb") as f:
        f.write(content.encode('utf-8'))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: \n")
        print("    python generate.py inlFilePath rtype nbits")
    else:
        rtype = sys.argv[2]
        nbits = int(sys.argv[3])
        if rtype != "int64_t":
            print("only import int64_t")
        generateConsts(sys.argv[1], rtype, nbits)
