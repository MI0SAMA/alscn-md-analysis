from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension("alscn_switch._xyzstream", ["src/alscn_switch/_xyzstream.c"]),
    ]
)
