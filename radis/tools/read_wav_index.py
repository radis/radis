import bisect
import json
import os

from radis.misc.utils import getProjectRoot

# Load index
wav_index_path = os.path.join(getProjectRoot(), "db", "wav_index.json")
with open(wav_index_path, "r") as f:
    wav_index = json.load(f)

wavno_keys = sorted(float(k) for k in wav_index.keys())

wavno_key_map = {float(k): k for k in wav_index.keys()}


def find_nearest_lower_wavno_index(wavno):
    """
    Find the index of the largest wav_no ≤ given wavno.
    If wavno is smaller than all, return index 0.
    """
    idx = bisect.bisect_right(wavno_keys, float(wavno))
    if idx == 0:
        return 0  # wavno smaller than any key, return smallest
    return idx - 1


def find_nearest_upper_wavno_index(wavno):
    """
    Find the index of the smallest wav_no ≥ given wavno.
    If wavno is larger than all, return last index.
    """
    idx = bisect.bisect_left(wavno_keys, float(wavno))
    if idx == len(wavno_keys):
        return len(wavno_keys) - 1  # wavno larger than any key, return largest
    return idx


def get_wavno_lower_offset(wavno):
    """
    Get the offset for the largest wav_no ≤ given wavno.
    Returns None if no such wavno exists.
    """
    index = find_nearest_lower_wavno_index(wavno)
    if index is not None:
        wavno_float = wavno_keys[index]
        wavno_key_str = wavno_key_map[wavno_float]
        print(
            f"Found lower wavno: {wavno_key_str} (float {wavno_float}) at index {index}"
        )
        return wav_index[wavno_key_str]
    else:
        print("No lower wavno found.")
        return None


def get_wavno_upper_offset(wavno):
    """
    Get the offset for the smallest wav_no ≥ given wavno.
    Returns None if no such wavno exists.
    """
    index = find_nearest_upper_wavno_index(wavno)
    if index is not None:
        wavno_float = wavno_keys[index]
        wavno_key_str = wavno_key_map[wavno_float]
        print(
            f"Found upper wavno: {wavno_key_str} (float {wavno_float}) at index {index}"
        )
        return wav_index[wavno_key_str]
    else:
        print("No upper wavno found.")
        return None


def offset_difference_from_lower_wavno(larger_wavno, lower_wavno):
    """
    Calculate the difference between two offsets:
    offset for smallest wavno ≥ larger_wavno and
    offset for largest wavno ≤ lower_wavno.

    Returns the absolute difference, or None if either offset is not found.
    """
    offset_upper = get_wavno_upper_offset(larger_wavno)
    offset_lower = get_wavno_lower_offset(lower_wavno)

    if offset_upper is None:
        print(f"No upper offset found for wavno: {larger_wavno}")
        return None
    if offset_lower is None:
        print(f"No lower offset found for wavno: {lower_wavno}")
        return None

    return abs(int(offset_upper) - int(offset_lower))
