#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_secstruct_design",
    version="0.1.0",
    description="Common secondary structure design algorithms for RNA",
    long_description=readme,
    long_description_content_type="text/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/rna_secstruct_design",
    packages=[
        "rna_secstruct_design",
    ],
    package_dir={"rna_secstruct_design": "rna_secstruct_design"},
    include_package_data=True,
    install_requires=requirements,
    python_requires=">=3.9",
    zip_safe=False,
    keywords="rna_secstruct_design",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    entry_points={
        "console_scripts": ["rna-struct-design = rna_secstruct_design.cli:cli"]
    },
)
