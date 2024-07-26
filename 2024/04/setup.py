from setuptools import setup
from Cython.Build import cythonize

setup(
    name='cython_arr',
    ext_modules=cythonize("cython_arr_nan.pyx"),
)
