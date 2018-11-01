# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py
import io
from setuptools import setup, find_packages


with io.open('README.md', "r", encoding="utf-8") as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='eth_radar',
    version='0.1.0',
    description='Echo top height estimation from radar data.',
    long_description=readme,
    author='Valentin Louf',
    author_email='valentin.louf@bom.gov.au',
    url='https://github.com/vlouf/eth_radar',
    license=license,
    packages=find_packages(exclude=('notebook'))
)
