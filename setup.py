#! /usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name="KinetiKit",
	version="0.1.8",
	author="Natalia Spitha",
	author_email="natalia.spitha@gmail.com",
    packages=find_packages(exclude=("tests", "tests.*")),
    python_requires=">=3.6",
    install_requires=[
        "matplotlib>=3.0",
        "numpy>=1.15.0",
        "scipy",
		"ipywidgets",
    ],
    description="Tools for comparing and kinetically simulating time-resolved data",
)
