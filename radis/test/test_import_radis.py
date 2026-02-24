# -*- coding: utf-8 -*-
"""Test radis import.

Start-up time can be monitored with Python 3.7:
https://dev.to/methane/how-to-speed-up-python-application-startup-time-nkf

See for instance:

https://stackoverflow.com/questions/16373510/improving-speed-of-python-module-import

which uses tuna

https://github.com/nschloe/tuna
"""

import json
import subprocess
import sys

import pytest

import radis

assert radis  # silences pyflakes


def _run_python(code):
    proc = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    return proc.stdout.strip()


@pytest.mark.fast
def test_import_radis_is_lazy():
    loaded = json.loads(
        _run_python(
            "import json, sys; import radis; "
            "print(json.dumps({"
            "'radis.misc.arrays': 'radis.misc.arrays' in sys.modules, "
            "'radis.tools.slit': 'radis.tools.slit' in sys.modules, "
            "'radis.lbl.broadening': 'radis.lbl.broadening' in sys.modules"
            "}))"
        )
    )
    assert loaded["radis.misc.arrays"] is False
    assert loaded["radis.tools.slit"] is False
    assert loaded["radis.lbl.broadening"] is False


@pytest.mark.fast
def test_import_radis_public_symbols():
    result = _run_python(
        "from radis import SpectrumFactory, calc_spectrum, config, load_spec; "
        "print(all(["
        "callable(calc_spectrum), "
        "isinstance(config, dict), "
        "callable(load_spec), "
        "isinstance(SpectrumFactory, type)"
        "]))"
    )
    assert result == "True"
