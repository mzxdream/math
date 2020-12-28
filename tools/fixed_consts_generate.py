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
        static constexpr {rtype} RAD2DEG = {rad2Deg}LL; //{rad2DegRaw}
        static constexpr {rtype} DEG2RAD = {deg2Rad}LL; //{deg2RadRaw}
        static constexpr {rtype} SQRT2 = {sqrt2}LL; //{sqrt2Raw}
        static constexpr {rtype} HALF_SQRT2 = {halfSqrt2}LL; //{halfSqrt2Raw}
        static constexpr {rtype} TWO_SQRT2 = {twoSqrt2}LL; //{twoSqrt2Raw}
        static constexpr {rtype} COMPARE_EPSILON = {compareEpsilon}LL; //{compareEpsilonRaw}
        static constexpr {rtype} ATAN2_P1 = {atan2p1}LL; //{atan2p1Raw}
        static constexpr {rtype} ATAN2_P2 = {atan2p2}LL; //{atan2p2Raw}
        static constexpr {rtype} ATAN2_P3 = {atan2p3}LL; //{atan2p3Raw}
        //2PI, 4PI, 8PI, 16PI, ...
        static constexpr {rtype} PI_TABLE[] = {{{piTable}}};
        //[0-90]
        static constexpr {rtype} COS_TABLE[] = {{{cosTable}}};
    }};
}} // namespace mzx

#endif"""


def generateConsts(inlFilePath, tbits, nbits):
    base = (1 << nbits)

    cosContent = ""
    cosTableCount = 3600
    for i in range(cosTableCount + 1):
        cosContent += ", {0}LL".format(
            round(math.cos(i * (math.pi / 2) / cosTableCount) * base))
    cosContent = cosContent[2:]

    piContent = ""
    i = 2
    while True:
        p = round(math.pi * i * base)
        if p > pow(2, tbits - 1) - 1:
            break
        piContent += ", {0}LL".format(p)
        i *= 2
    piContent = piContent[2:]

    compareEpsilonRaw = pow(0.1, int(nbits/4)) * 100
    rtype = "int{}_t".format(tbits)
    datas = {
        "rtype": rtype,
        "nbits": nbits,
        "pi": round(math.pi * base),
        "piRaw": math.pi,
        "halfPI": round(math.pi * base / 2.0),
        "halfPIRaw": math.pi / 2.0,
        "twoPI": round(math.pi * base * 2.0),
        "twoPIRaw": math.pi * 2.0,
        "rad2Deg": round(180 * base / math.pi),
        "rad2DegRaw": 180 / math.pi,
        "deg2Rad": round(math.pi * base / 180),
        "deg2RadRaw": math.pi / 180,
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


def main():
    if len(sys.argv) < 3:
        print("usage: \n")
        print("    python generate.py inlFilePath rtype nbits")
    else:
        tbits = int(sys.argv[2])
        nbits = int(sys.argv[3])
        if tbits != 16 and tbits != 32 and tbits != 64:
            print("only import 16/32/64")
            return
        if nbits < 1 or nbits > tbits - 10:# 360deg limit
            print("nbits need >= 1 and <= tbits - 10")
            return
        generateConsts(sys.argv[1], tbits, nbits)

if __name__ == "__main__":
    main()
