import json
import os
from bisect import bisect_left, bisect_right

from radis.misc.utils import getProjectRoot


def key_pairs(load_min, load_max, key_pos=0, wrap=False):
    # extract numeric keys
    wav_index_path = os.path.join(getProjectRoot(), "db", "wav_index.json")
    with open(wav_index_path, "r") as f:
        raw = json.load(f)

    # raw["data"] maps string keys to numeric offsets
    wav_index = {float(k): int(v) for k, v in raw["data"].items()}

    if isinstance(wav_index, dict):
        keys = [float(k) for k in wav_index.keys()]
    else:
        if not wav_index:
            return []
        first = wav_index[0]
        if isinstance(first, (list, tuple)):
            keys = [item[key_pos] for item in wav_index]
            # convert string keys to float if necessary
            keys = [float(k) if isinstance(k, str) else k for k in keys]
        else:
            keys = [float(k) if isinstance(k, str) else k for k in wav_index]

    keys = sorted(keys)
    if not keys:
        return []

    s = bisect_right(keys, load_min) - 1
    if s < 0:
        s = 0
    e = bisect_left(keys, load_max)
    if e >= len(keys):
        e = len(keys) - 1

    pairs = []
    if s < e:
        for i in range(s, e):
            pairs.append((keys[i], keys[i + 1]))

    if wrap and keys:
        pairs.append((keys[e], keys[s]))

    return pairs
