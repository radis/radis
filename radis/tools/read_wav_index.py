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


def get_wavno_for_byte_range(start_byte: int, end_byte: int) -> tuple[float, float]:
    """
    Return (start_wavno, end_wavno) for the given byte range.
    Requires: wavno_keys (sorted list of float keys) and wav_index (dict mapping key->offset).
    """
    if start_byte > end_byte:
        start_byte, end_byte = end_byte, start_byte

    wav_offsets = [wav_index[k] for k in wavno_keys]

    # largest offset <= start_byte
    i = bisect.bisect_right(wav_offsets, start_byte)
    start_idx = 0 if i == 0 else i - 1

    # smallest offset >= end_byte
    j = bisect.bisect_left(wav_offsets, end_byte)
    end_idx = len(wav_offsets) - 1 if j == len(wav_offsets) else j

    return wavno_keys[start_idx], wavno_keys[end_idx]
