#! /usr/bin/env python3

from setuptools import setup

setup(
    name="KinetiKit",
    packages=['KinetiKit'],
    python_requires=">=3.6",
    install_requires=[
        "matplotlib>=3.0",
        "numpy>=1.15.0",
        "scipy",
    ],
    version='0.1.0',
    description="Tools for comparing and kinetically simulating time-resolved data",
)
