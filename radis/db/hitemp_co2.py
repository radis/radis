"""
HITEMP CO2 Block-Aligned Partial Decompression.

Implements efficient partial downloading and decompression of the HITEMP CO2
database using block-aligned bzip2 decompression for specific wavenumber ranges.

@author: dcmvd
"""

import bz2
import io
import os

import numpy as np
from tqdm import tqdm

from radis.misc.utils import getProjectRoot

# Use absolute path to avoid issues
project_root = getProjectRoot()
offset_path = os.path.join(project_root, "db", "offset_arr.npy")

offsets, pads = np.load(offset_path)
sizes = offsets[1:] - offsets[:-1]


pad_bufs = {
    65: b"BZh91AY&SY Q\xb2\x1c\x00\x00\x03Z\x00\x00\x10\x0c\x00@\x00\x00\n \x000\xc0\x084\xf2 b\xfd",
    130: b"BZh91AY&SY\x89\xbf\xa3\xf5\x00\x00\x03Z\x00\x00\x10\x0e\x00 \x00\x00\n \x001\x0c\x01\x06\x99\xa1?\x19\n",
    5: b"BZh91AY&SYY\x9e\xb7A\x00\x00\x03Z\x00\x00\x10\x0c\x00\x10\x00\x00\n \x000\xc0\x08a\xa1gPx",
    10: b"BZh91AY&SY\x13\xc2\x98\xc5\x00\x00\x03Z\x00\x00\x10\x0c\x00\x08\x00\x00\n \x000\xc0\x08a\xa1e@\xf1",
    21: b'BZh91AY&SY\xe7%cn\x00\x00\x03\xda\x00@\x10\x08\x00\x04\x00\x00\n \x00"\x18h0\x06/\xa0c',
    43: b"BZh91AY&SYm\x96(6\x00\x00\x03Z\x00\x00\x10\x0e\x00\x02\x00\x00\n \x001\x0c\x01\x01\xb2\x89\xc6#\xe6",
    86: b'BZh91AY&SY|\x80\xf4\xd9\x00\x00\x03Z\x00\x00\x10\x0c\x00\x01\x00\x00\n \x00"\x18h0\x02\xcf)\x8c',
    172: b"BZh91AY&SY\xc9\xec\x1b[\x00\x00\x03Z\x00\x00\x10\x0e\x00\x00\x80\x00\n \x001\x0c\x01\r1\xa8Y0\xa7\x98",
}


def get_bz2(session, file_url, offset=None, size=None, verbose=True):
    """
    Download a byte range from a remote bzip2 file using HTTP range requests.

    Parameters
    ----------
    session : requests.Session
        Authenticated session
    file_url : str
        URL of the remote file
    offset, size : int, optional
        Byte range to download. If None, downloads entire file.
    verbose : bool, default True
        If True, prints progress and status messages

    Returns
    -------
    bytes
        Downloaded file content
    """
    if offset is not None and size is not None:
        range_headers = {"Range": f"bytes={offset}-{offset + size - 1}"}
    else:
        range_headers = {}

    # Download the file chunk using the authenticated session
    with session.get(file_url, headers=range_headers, stream=True) as r:
        total = int(r.headers.get("content-length", 0))

        buf = io.BytesIO()
        with tqdm(
            desc="Downloading",
            total=total,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            disable=not verbose,
        ) as bar:
            for chunk in r.iter_content(chunk_size=8192):
                buf.write(chunk)
                bar.update(len(chunk))

    if verbose:
        print("Chunk Download complete")
    return buf.getvalue()


def partial_download_co2_chunk(
    target_wn_min, target_wn_max, session, output_file_path, verbose=True
):
    """
    Download and decompress a specific wavenumber range from HITEMP CO2 database.

    Uses block-aligned partial decompression to extract only the requested wavenumber range without downloading the entire ~6GB file.

    Parameters
    ----------
    target_wn_min, target_wn_max : float
        Wavenumber range to extract (cm⁻¹)
    session : requests.Session
        Authenticated HITRAN session
    output_file_path : str
        Where to save the decompressed data
    verbose : bool, default True
        If True, prints progress and status messages
    """

    wavenumber_path = os.path.join(project_root, "db", "wavenumber_arr.npy")
    wavenumbers = np.load(wavenumber_path)
    i_min = np.searchsorted(wavenumbers, target_wn_min, side="left")
    i_max = np.searchsorted(wavenumbers, target_wn_max, side="right")

    # Include adjacent blocks to ensure complete coverage
    i_min = max(0, min(i_min - 1, len(offsets) - 1))
    i_max = min(i_max + 1, len(offsets) - 1)

    if verbose:
        print(f"Found target in compressed blocks: {i_min} to {i_max}")
        print(
            f"Wavenumber range: {wavenumbers[i_min]:.6f} to {wavenumbers[i_max]:.6f} cm-1"
        )

    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    offset = offsets[i_min]
    size = offsets[i_max + 1] - offsets[i_min]

    file_url = "https://hitran.org/files/HITEMP/bzip2format/02_HITEMP2024.par.bz2"
    buf = get_bz2(session, file_url, offset=offset, size=size, verbose=verbose)

    # Perform block-aligned bzip2 decompression
    pad_index_to_key = [65, 130, 5, 10, 21, 43, 86, 172]

    decomp = bz2.BZ2Decompressor()
    pad_key = pad_index_to_key[pads[i_min]]
    decomp.decompress(pad_bufs[pad_key])
    raw = decomp.decompress(buf)
    decomp = None
    buf = None

    # Trim partial lines to ensure valid HITRAN format
    first_nl = raw.find(b"\n")
    last_nl = raw.rfind(b"\n")
    if first_nl != -1 and last_nl != -1 and first_nl < last_nl:
        data = raw[first_nl + 1 : last_nl]
    else:
        data = raw

    with open(output_file_path, "wb") as fw:
        fw.write(data)
    if verbose:
        print(f"Written decompressed data to {output_file_path}")
