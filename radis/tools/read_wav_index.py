import bisect
import json
import os
import pickle

import indexed_bzip2 as ibz2

from radis.misc.utils import getProjectRoot

# Load index
wav_index_path = os.path.join(getProjectRoot(), "db", "wav_index.json")
with open(wav_index_path, "r") as f:
    wav_index = json.load(f)

wavno_keys = sorted(float(k) for k in wav_index.keys())

wavno_key_map = {float(k): k for k in wav_index.keys()}


def find_nearest_lower_wavno_index(wavno):
    """
    Find the index of largest wav_no ≤ given wavno using binary search.
    Returns None if no such wavno exists.
    """
    idx = bisect.bisect_right(wavno_keys, float(wavno))
    return idx - 1 if idx > 0 else None


def find_nearest_upper_wavno_index(wavno):
    """
    Find the index of smallest wav_no ≥ given wavno using binary search.
    Returns None if no such wavno exists.
    """
    idx = bisect.bisect_left(wavno_keys, float(wavno))
    return idx if idx < len(wavno_keys) else None


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


def read_and_write_chunked_for_CO2(
    bz2_file_path,
    index_file_path,
    start_offset,
    bytes_to_read,
    chunk_size=500 * 1024 * 1024,
    output_prefix="CO2_HITEMP",
):
    """
    Read `bytes_to_read` bytes from a bzip2 file starting at `start_offset`, in chunks, and
    write each chunk to its own file named:
        {output_prefix}_{start_mb}mb.par
    where `start_mb` is the starting offset of that chunk in megabytes.

    """
    # Load block offsets
    with open(index_file_path, "rb") as f:
        block_offsets = pickle.load(f)

    # Open and prepare the bzip2 file
    f = ibz2.open(bz2_file_path, parallelization=os.cpu_count())
    f.set_block_offsets(block_offsets)
    f.seek(start_offset)

    total_read = 0
    # Read in chunks and write separate files
    while total_read < bytes_to_read:
        to_read = min(chunk_size, bytes_to_read - total_read)
        data = f.read(to_read)
        if not data:
            break  # end of file

        # Compute current offset for naming
        current_offset = start_offset + total_read
        start_mb = current_offset // (1024 * 1024)
        out_name = f"{output_prefix}_{start_mb}mb.par"

        # Write this chunk
        with open(out_name, "wb") as out_file:
            out_file.write(data)

        total_read += len(data)
        print(
            f"Wrote chunk {start_mb}MiB → {len(data)} bytes to {out_name} ({total_read}/{bytes_to_read})"
        )

    f.close()
    print(
        f"Finished: wrote {total_read} bytes split into { (total_read + chunk_size - 1) // chunk_size } file(s). "
    )


read_and_write_chunked_for_CO2(
    "02_HITEMP2024.par.bz2",
    "CO2_indexed_offsets.dat",
    500 * 1024 * 1024,
    2000 * 1024 * 1024,
)
