import bisect
import json
import os

from radis.misc.utils import getProjectRoot

# Load index and convert keys to float
wav_index_path = os.path.join(getProjectRoot(), "db", "wav_index.json")
with open(wav_index_path, "r") as f:
    raw = json.load(f)

# raw["data"] maps string keys to numeric offsets
wav_index = {float(k): int(v) for k, v in raw["data"].items()}
wavno_keys = sorted(wav_index)


def find_nearest_lower_wavno_index(wavno: float) -> int:
    """
    Find the index of the largest wav_no ≤ given wavno.
    If wavno is smaller than all, return index 0.
    """
    idx = bisect.bisect_right(wavno_keys, wavno)
    if idx == 0:
        return 0  # wavno smaller than any key, return smallest
    return idx - 1


def find_nearest_upper_wavno_index(wavno: float) -> int:
    """
    Find the index of the smallest wav_no ≥ given wavno.
    If wavno is larger than all, return last index.
    """
    idx = bisect.bisect_left(wavno_keys, wavno)
    if idx == len(wavno_keys):
        return len(wavno_keys) - 1  # wavno larger than any key, return largest
    return idx


def get_wavno_lower_offset(wavno: float) -> int | None:
    """
    Get the offset for the largest wav_no ≤ given wavno.
    Returns None if no such wavno exists.
    """
    index = find_nearest_lower_wavno_index(wavno)
    key = wavno_keys[index]
    return wav_index.get(key)


def get_wavno_upper_offset(wavno: float) -> int | None:
    """
    Get the offset for the smallest wav_no ≥ given wavno.
    Returns None if no such wavno exists.
    """
    index = find_nearest_upper_wavno_index(wavno)
    key = wavno_keys[index]
    return wav_index.get(key)


def offset_difference_from_lower_wavno(
    larger_wavno: float, lower_wavno: float
) -> int | None:
    """
    Calculate the difference between two offsets:
    offset for smallest wavno ≥ larger_wavno and
    offset for largest wavno ≤ lower_wavno.

    Returns the absolute difference, or None if either offset is not found.
    """
    upper = get_wavno_upper_offset(larger_wavno)
    lower = get_wavno_lower_offset(lower_wavno)

    if upper is None:
        return None
    if lower is None:
        return None

    return abs(upper - lower)
