# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py
import io
from setuptools import setup, find_packages


with io.open('README.md', "r", encoding="utf-8") as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

# What packages are required for this module to be executed?
REQUIRED = [
    "numpy", "numba", "scipy", "arm_pyart"
]

setup(
    name='echotop',
    version='1.0',
    description='Echo top height estimation from radar data.',
    long_description=readme,
    author='Valentin Louf',
    author_email='valentin.louf@monash.edu',
    url='https://github.com/vlouf/eth_radar',
    license=license,
    install_requires=REQUIRED,
    include_package_data=True,
    packages=find_packages(exclude=('notebook')),
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',        
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Intended Audience :: Science/Research',
    ],
)
