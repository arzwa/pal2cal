from setuptools import setup

import pal2cal

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pal2cal",
    version=pal2cal.__version__,
    author=pal2cal.__author__,
    author_email="arzwa@psb.vib-ugent.be",
    py_modules=["pal2cal"],
    description="Yet another pal2nal script",
    license="MIT",
    python_requires='>=3.5',
    entry_points='''
        [console_scripts]
        pal2cal=pal2cal:main
        ''',
    install_requires=[
        'click>=7.0',
        'biopython>=1.75',])
