from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext
from numpy import __path__ as npypath

ext_modules = [
    Extension(
        "cython_radis",
        sources=["cython_radis.pyx"],
        include_dirs=[npypath[0] + "/core/include"],
        language="c",
        extra_compile_args=["/O2", "/favor:INTEL64", "/fp:fast"],
        extra_link_args=[],
    )
]
setup(name="cython_radis", cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)
