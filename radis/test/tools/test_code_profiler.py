# -*- coding: utf-8 -*-
"""
Test the :py:class:`~radis.tools.code_profiler.CodeProfiler` class.
"""

import os
import sys
import tempfile

import pytest


@pytest.mark.fast
def test_code_profiler_init():
    """Test CodeProfiler initialization and dump."""
    from radis.tools.code_profiler import CodeProfiler

    original_profile = sys.getprofile()
    temp_file = None

    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            temp_file = f.name

        profiler = CodeProfiler(filename=temp_file)

        assert profiler.white_list == {"radis"}
        assert profiler.exclusions == {"<"}

        profiler.dump(filename=temp_file)
        assert os.path.exists(temp_file)

    finally:
        sys.setprofile(original_profile)
        if temp_file and os.path.exists(temp_file):
            os.remove(temp_file)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
