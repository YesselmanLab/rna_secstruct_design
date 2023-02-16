#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()


with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name='rna_secstruct_design',
    version='0.1.0',
    description='A set of scripts for handdling common secondary structure design algorithms',
    long_description=readme,
    long_description_content_type="test/markdown",
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/rna_secstruct_design',
    packages=[
        'rna_secstruct_design',
    ],
    package_dir={'rna_secstruct_design': 'rna_secstruct_design'},
    py_modules=[
        'rna_secstruct_design/cli'
        'rna_secstruct_design/constraints',
        'rna_secstruct_design/helix_randomizer',
        'rna_secstruct_design/logger',
        'rna_secstruct_design/mutations',
        'rna_secstruct_design/selection',
        'rna_secstruct_design/util',
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='rna_secstruct_design',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
            'rna-struct-design = rna_secstruct_design.cli:cli'
        ]
    }
)
