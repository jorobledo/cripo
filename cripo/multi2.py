"""Bragg-cut generator translated from MULTI2."""

from __future__ import annotations

import math
from dataclasses import dataclass

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


def _unpack_two_digit_table(rows: list[int]) -> list[int]:
    """Expand legacy packed two-digit table rows into a flat Python list."""
    values: list[int] = []
    for row in rows:
        digits = f"{row:d}"
        if len(digits) % 2 == 1:
            digits = "0" + digits
        values.extend(int(digits[i : i + 2]) for i in range(0, len(digits), 2))
    return values


ICUB_VALUES = _unpack_two_digit_table(ICUB)
I1_VALUES = _unpack_two_digit_table(I1)
IWHEX_VALUES = _unpack_two_digit_table(IWHEX)
ISTETR_VALUES = _unpack_two_digit_table(ISTETR)
IWTETR_VALUES = _unpack_two_digit_table(IWTETR)

HEXAGONAL_TYPES = {1, 2, 6}
CUBIC_TYPES = {3, 4, 5}
GRAPHITE_TYPE = 1
HEXAGONAL_CLOSE_PACKED_TYPE = 2
HEXAGONAL_SPECIAL_TYPE = 6
BODY_CENTERED_CUBIC_TYPE = 3
FACE_CENTERED_CUBIC_TYPE = 4
TETRAGONAL_TYPE = 7
BRAGG_TABLE_OVERFLOW = 103


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
    if indic in HEXAGONAL_TYPES:
        return _multi2_hex(a, c, elim, indic, idim, ier)
    if indic in CUBIC_TYPES:
        return _multi2_cub(a, elim, indic, idim, ier)
    if indic == TETRAGONAL_TYPE:
        return _multi2_tetr(a, c, elim, idim, ier)
    raise ValueError(f"Unsupported crystal type: {indic}")


def _sort_pairs(a1: list[float], b1: list[float]) -> tuple[list[float], list[float]]:
    """Sort paired energy and weight arrays by ascending energy."""
    pairs = sorted(zip(a1, b1), key=lambda item: item[0])
    return [p[0] for p in pairs], [p[1] for p in pairs]


def _find_table_limit(values: list[int], maximum: int, default_limit: int, fallback_limit: float) -> tuple[int, float]:
    """Return the usable table length and the capped search range."""
    for index, value in enumerate(values[:default_limit], start=1):
        if value > maximum:
            return index - 1, float(maximum)
    return default_limit, fallback_limit


def _append_peak(energies: list[float], weights: list[float], energy: float, weight: float) -> None:
    """Append one Bragg-cut peak to the output lists."""
    energies.append(energy)
    weights.append(weight)


def _hexagonal_atom_divisor(crystal_type: int) -> float:
    """Return the multiplicity divisor for the selected hexagonal structure."""
    if crystal_type == HEXAGONAL_CLOSE_PACKED_TYPE:
        return 2.0
    if crystal_type == HEXAGONAL_SPECIAL_TYPE:
        return 1.0
    return 4.0


def _hexagonal_l3_rules(crystal_type: int, base_index: int) -> tuple[int, int, int]:
    """Return stepping and parity factors for the hexagonal/turbostratic loops."""
    if crystal_type == HEXAGONAL_SPECIAL_TYPE:
        return 1, 1, 1

    parity_factor = 16 if crystal_type == GRAPHITE_TYPE else 4
    if base_index % 3 == 0:
        return 2, parity_factor, 3
    return 1, 1, 3


def _cubic_weight_scale(crystal_type: int) -> float:
    """Return the crystal-type-dependent multiplicity scale for cubic lattices."""
    if crystal_type == BODY_CENTERED_CUBIC_TYPE:
        return 2.0
    if crystal_type == FACE_CENTERED_CUBIC_TYPE:
        return 4.0
    return 1.0


def _allowed_cubic_index(is_value: int, crystal_type: int) -> bool:
    """Check whether a cubic reflection index is allowed for the structure."""
    if crystal_type == 5:
        return True
    if crystal_type == BODY_CENTERED_CUBIC_TYPE:
        return is_value % 2 == 0
    return is_value % 8 in (0, 3, 4)


def _tetragonal_component(value: int) -> int:
    """Classify an index by the mod-4 parity bucket used in the original code."""
    if value % 2 != 0:
        return 1
    if value % 4 != 0:
        return 2
    return 3


def _multi2_hex(a: float, c: float, elim: float, indic: int, idim: int, ier: int) -> Multi2Result:
    """Handle hexagonal and graphite-like crystal structures."""
    energy_scale = (0.286015e-8 / a) ** 2 / 3.0
    scaled_limit = elim / energy_scale
    atom_divisor = _hexagonal_atom_divisor(indic)
    table_limit, scaled_limit = _find_table_limit(I1_VALUES, int(scaled_limit), default_limit=36, fallback_limit=97.0)

    axial_ratio = 3.0 * a * a / 4.0 / c / c
    v0 = c * a * a * math.sqrt(3.0) / 2.0
    energies: list[float] = []
    weights: list[float] = []

    for base_index, base_weight in zip(I1_VALUES[:table_limit], IWHEX_VALUES[:table_limit]):
        if scaled_limit < base_index:
            continue
        l3_step, even_factor, odd_factor = _hexagonal_l3_rules(indic, base_index)
        l3_limit = int(math.sqrt((scaled_limit - base_index) / axial_ratio) + 1.0)
        for l3 in range(0, l3_limit, l3_step):
            scaled_energy = base_index + axial_ratio * l3 * l3
            if scaled_energy <= 0.0 or scaled_energy > scaled_limit:
                continue
            structure_factor = even_factor if l3 % 2 == 0 else odd_factor
            multi = base_weight * structure_factor
            if l3 == 0:
                multi //= 2
            _append_peak(energies, weights, scaled_energy * energy_scale, multi / atom_divisor)
            if len(energies) >= idim:
                return Multi2Result(energies, weights, energies[-1], len(energies), BRAGG_TABLE_OVERFLOW, v0)

    energies, weights = _sort_pairs(energies, weights)
    last = energies[-1] if energies else 0.0
    return Multi2Result(energies, weights, last, len(energies), ier, v0)


def _multi2_cub(a: float, elim: float, indic: int, idim: int, ier: int) -> Multi2Result:
    """Handle simple cubic, BCC, and FCC crystal structures."""
    v0 = a * a * a
    energy_scale = (0.286015e-8 / a / 2.0) ** 2
    max_index = min(int(elim / energy_scale), 800)
    weight_scale = _cubic_weight_scale(indic)

    energies: list[float] = []
    weights: list[float] = []
    for is_value in range(1, max_index + 1):
        if not _allowed_cubic_index(is_value, indic):
            continue
        multiplicity_index = ICUB_VALUES[is_value - 1]
        if multiplicity_index >= len(MCUB) or MCUB[multiplicity_index] == 0:
            continue
        _append_peak(energies, weights, is_value * energy_scale, MCUB[multiplicity_index] * weight_scale)
        if len(energies) >= idim:
            break

    last = energies[-1] if energies else 0.0
    return Multi2Result(energies, weights, last, len(energies), ier, v0)


def _multi2_tetr(a: float, c: float, elim: float, idim: int, ier: int) -> Multi2Result:
    """Handle the tetragonal tin case from the original implementation."""
    energy_scale = (0.5 * 0.286015e-8 / a) ** 2
    scaled_limit = elim / energy_scale
    table_limit, scaled_limit = _find_table_limit(ISTETR_VALUES, int(scaled_limit), default_limit=43, fallback_limit=98.0)

    axial_ratio = a * a / c / c
    v0 = a * a * c
    energies: list[float] = []
    weights: list[float] = []

    for base_index, base_weight in zip(ISTETR_VALUES[:table_limit], IWTETR_VALUES[:table_limit]):
        if scaled_limit < base_index:
            continue
        base_component = _tetragonal_component(base_index)
        l3_limit = int(math.sqrt((scaled_limit - base_index) / axial_ratio) + 1.0)
        for l3 in range(l3_limit):
            l3_component = _tetragonal_component(l3)
            if l3_component != base_component:
                continue
            iu = 4 if l3_component == 1 else 8
            if l3 == 0:
                iu //= 2
            scaled_energy = base_index + axial_ratio * l3 * l3
            if scaled_energy <= 0.0 or scaled_energy > scaled_limit:
                continue
            _append_peak(energies, weights, scaled_energy * energy_scale, base_weight * iu)
            if len(energies) >= idim:
                return Multi2Result(energies, weights, energies[-1], len(energies), BRAGG_TABLE_OVERFLOW, v0)

    energies, weights = _sort_pairs(energies, weights)
    last = energies[-1] if energies else 0.0
    return Multi2Result(energies, weights, last, len(energies), ier, v0)
