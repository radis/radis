import json
import os
from bisect import bisect_left, bisect_right
from typing import List, Tuple

from radis.misc.utils import getProjectRoot


def get_key_pairs(
    min_wavenumber: float, max_wavenumber: float
) -> List[Tuple[float, float]]:
    """
    Read db/wav_index.json and return consecutive (wavenumber_i, wavenumber_{i+1}) pairs
    whose intervals overlap the range [min_wavenumber, max_wavenumber].
    """
    wav_index_path = os.path.join(getProjectRoot(), "db", "wav_index.json")
    if not os.path.exists(wav_index_path):
        raise FileNotFoundError(f"{wav_index_path} not found")

    with open(wav_index_path, "r") as f:
        raw = json.load(f)

    data_dict = raw["data"]

    wavenumbers = sorted(float(k) for k in data_dict.keys())

    if not wavenumbers:
        return []

    # Find index range that overlaps [min_wavenumber, max_wavenumber]
    start_idx = bisect_right(wavenumbers, min_wavenumber) - 1
    if start_idx < 0:
        start_idx = 0

    end_idx = bisect_left(wavenumbers, max_wavenumber)
    if end_idx >= len(wavenumbers):
        end_idx = len(wavenumbers) - 1

    pairs: List[Tuple[float, float]] = []
    if start_idx < end_idx:
        for i in range(start_idx, end_idx):
            pairs.append((wavenumbers[i], wavenumbers[i + 1]))

    return pairs
