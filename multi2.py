"""Bragg-cut generator translated from MULTI2."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .cs1 import cs1

ICUB = [
    2030402, 4010505, 3050506, 2010905, 5050809, 5010509, 1071206, 4010912, 6090909, 5011205, 5050915, 3010912,
    5091410, 9011512, 1120509, 2011512, 9051515, 8010915, 5111909, 5010915, 9120916, 5011909, 1091922, 5011509,
    6121809, 12010924, 7121215, 1012212, 12091515, 9011219, 1091913, 4012222, 9190922, 9012409, 9121515, 6010915,
    5112615, 12011924, 1151522, 5011512, 15051826, 9011915, 5192213, 5012224, 12191509, 9012219, 1092226, 3010926,
    10153109, 14011919, 9151222, 5013115, 12121531, 15011215, 1151915, 9012226, 9191523, 5013119, 12092822, 1011522,
    12172122, 15012224, 1242222, 2012609, 15151526, 15011228, 5093115, 9012633, 15192215, 9012412, 1121531, 8012226,
    9153116, 19012428, 11261226, 5012215, 15123126, 9011522, 1223615, 5012219, 16151531, 9013324, 12121535, 9012219,
    9222515, 19012838, 1191915, 5013124, 22152626, 19011215, 9172820, 1012231, 9262235, 15013819, 1122626, 5011531,
    9193115, 18012633, 12261922, 6013112, 24092631, 9012628, 1222626, 12013131, 15281522, 12011919, 12093538, 7012226,
    12153522, 22011531, 1191534, 1014015, 15223231, 15012422, 9263615, 12011539, 19241526, 12013126, 1193122, 9011528,
    13153126, 19013631, 9282431, 1013822, 22121535, 22012426, 1173815, 4012639, 22222626, 9012622, 19191940, 9012231,
    9283823, 24013122, 1281922, 9013519, 15123540, 15012439, 12153315, 9013528, 15281526, 9013922, 1192640, 6012626,
    15262715, 26012638, 11242835, 5013126, 24154035, 19011522, 1193531, 12012236, 22312226, 15014119, 15192222, 1012240,
    12154223, 15013142, 1312431, 5013124, 26122638, 18012626, 5263931, 15013824, 15331531, 19013328, 1153935, 9012222,
    13293522, 22013139, 19192231, 5014022, 24192643, 22011935, 1313615, 5012641, 9223535, 15014019, 19153535, 12012226,
    19153026, 22013339, 1352237, 9014324, 26193131, 22011535, 9283626, 1013540, 26353122, 9013319, 1152243, 3013531,
    9224031, 31013936, 15312426, 10013531, 19192626, 19012640, 1154415, 14012639,
]
MCUB = [
    0, 6, 8, 12, 24, 30, 32, 36, 48, 54, 56, 72,
    78, 84, 96, 102, 104, 108, 120, 126, 132, 144, 150, 168,
    180, 192, 204, 216, 224, 228, 240, 252, 264, 270, 288, 312,
    318, 336, 360, 384, 408, 432, 480, 504,
]
I1 = [4030100, 13120907, 25211916, 36312827, 48433937, 61575249, 73676463, 81797675, 97939184]
IWHEX = [12121202, 24121224, 12242412, 12242412, 12242424, 24242436, 24241224, 12242412, 24244824]
ISTETR = [4020100, 10090805, 18171613, 29262520, 37363432, 49454140, 58535250, 68656461, 80747372, 89858281, 989790]
IWTETR = [4040401, 8040408, 4080408, 8081208, 8040804, 4080808, 8080812, 8160408, 8080804, 8160804, 40808]


@dataclass
class Multi2Result:
    """Container for the Bragg-cut energies and multiplicity factors."""

    a1: list[float]
    b1: list[float]
    elim: float
    idim: int
    ier: int
    v0: float


def multi2(a: float, c: float, elim: float, indic: int, idim: int) -> Multi2Result:
    """Generate Bragg-cut energies and weights for the requested crystal type."""
    ier = 0
    if indic in (1, 2, 6):
        return _multi2_hex(a, c, elim, indic, idim, ier)
    if indic in (3, 4, 5):
        return _multi2_cub(a, elim, indic, idim, ier)
    if indic == 7:
        return _multi2_tetr(a, c, elim, idim, ier)
    raise ValueError(f"Unsupported crystal type: {indic}")


def _sort_pairs(a1: list[float], b1: list[float]) -> tuple[list[float], list[float]]:
    """Sort paired energy and weight arrays by ascending energy."""
    pairs = sorted(zip(a1, b1), key=lambda item: item[0])
    return [p[0] for p in pairs], [p[1] for p in pairs]


def _multi2_hex(a: float, c: float, elim: float, indic: int, idim: int, ier: int) -> Multi2Result:
    """Handle hexagonal and graphite-like crystal structures."""
    fact = 0.286015e-8 / a
    fact = fact * fact / 3.0
    f = elim / fact
    an = 4.0
    if indic == 2:
        an = 2.0
    if indic == 6:
        an = 1.0

    ismax = int(f)
    islim = 36
    for i in range(1, 37):
        j11, _ = cs1(I1, i, 0)
        if j11 > ismax:
            islim = i - 1
            break
    else:
        f = 97.0

    ak = 3.0 * a * a / 4.0 / c / c
    v0 = c * a * a * math.sqrt(3.0) / 2.0
    out_a1: list[float] = []
    out_b1: list[float] = []

    for iis in range(1, islim + 1):
        is_value, _ = cs1(I1, iis, 0)
        iww, _ = cs1(IWHEX, iis, 0)
        if f < is_value:
            continue
        l3lim = int(math.sqrt((f - is_value) / ak) + 1.0)
        iu = 1
        l3par = 1
        l3non = 1
        if indic != 6:
            l3non = 3
            is1 = is_value - (is_value // 3) * 3
            if is1 == 0:
                iu = 2
                l3par = 16 if indic == 1 else 4
        for ll3 in range(1, l3lim + 1, iu):
            l3 = ll3 - 1
            energy = is_value + ak * l3 * l3
            if energy <= 0.0 or energy > f:
                continue
            iform = l3par if l3 % 2 == 0 else l3non
            multi = iww * iform
            if l3 == 0:
                multi //= 2
            out_a1.append(energy * fact)
            out_b1.append(multi / an)
            if len(out_a1) >= idim:
                return Multi2Result(out_a1, out_b1, out_a1[-1], len(out_a1), 103, v0)

    out_a1, out_b1 = _sort_pairs(out_a1, out_b1)
    last = out_a1[-1] if out_a1 else 0.0
    return Multi2Result(out_a1, out_b1, last, len(out_a1), ier, v0)


def _multi2_cub(a: float, elim: float, indic: int, idim: int, ier: int) -> Multi2Result:
    """Handle simple cubic, BCC, and FCC crystal structures."""
    v0 = a * a * a
    fact = (0.286015e-8 / a / 2.0) ** 2
    islim = min(int(elim / fact), 800)
    f2dn = 1.0
    if indic == 3:
        f2dn = 2.0
    if indic == 4:
        f2dn = 4.0

    out_a1: list[float] = []
    out_b1: list[float] = []
    for is_value in range(1, islim + 1):
        if indic != 5:
            if indic == 3 and is_value % 2 != 0:
                continue
            if indic != 3:
                iresto = is_value % 8
                if iresto not in (0, 3, 4):
                    continue
        j1, _ = cs1(ICUB, is_value, 0)
        if j1 >= len(MCUB) or MCUB[j1] == 0:
            continue
        out_a1.append(is_value * fact)
        out_b1.append(MCUB[j1] * f2dn)
        if len(out_a1) >= idim:
            break

    last = out_a1[-1] if out_a1 else 0.0
    return Multi2Result(out_a1, out_b1, last, len(out_a1), ier, v0)


def _multi2_tetr(a: float, c: float, elim: float, idim: int, ier: int) -> Multi2Result:
    """Handle the tetragonal tin case from the original implementation."""
    fact = (0.5 * 0.286015e-8 / a) ** 2
    f = elim / fact
    ismax = int(f)
    islim = 43
    for i in range(1, 44):
        j11, _ = cs1(ISTETR, i, 0)
        if j11 > ismax:
            islim = i - 1
            break
    else:
        f = 98.0

    ak = a * a / c / c
    v0 = a * a * c
    out_a1: list[float] = []
    out_b1: list[float] = []

    for iis in range(1, islim + 1):
        is_value, _ = cs1(ISTETR, iis, 0)
        iww, _ = cs1(IWTETR, iis, 0)
        icomp1 = 1
        if is_value % 2 == 0:
            icomp1 = 2
            if is_value % 4 == 0:
                icomp1 = 3
        if f < is_value:
            continue
        l3lim = int(math.sqrt((f - is_value) / ak) + 1.0)
        for ll3 in range(1, l3lim + 1):
            l3 = ll3 - 1
            icomp2 = 1
            if l3 % 2 == 0:
                icomp2 = 2
                if l3 % 4 == 0:
                    icomp2 = 3
            if icomp2 != icomp1:
                continue
            iu = 4 if icomp2 == 1 else 8
            if l3 == 0:
                iu //= 2
            energy = is_value + ak * l3 * l3
            if energy <= 0.0 or energy > f:
                continue
            out_a1.append(energy * fact)
            out_b1.append(iww * iu)
            if len(out_a1) >= idim:
                return Multi2Result(out_a1, out_b1, out_a1[-1], len(out_a1), 103, v0)

    out_a1, out_b1 = _sort_pairs(out_a1, out_b1)
    last = out_a1[-1] if out_a1 else 0.0
    return Multi2Result(out_a1, out_b1, last, len(out_a1), ier, v0)
