#! /usr/bin/env python
from setuptools import setup, find_packages
from oceanwaves import __version__


setup(
    name='oceanwaves',
    version=__version__,
    author='Patricia Wiberg, Christopher Sherwood',
    author_email='pw3c@virginia.edu',
    url='http://github.com/csdms/oceanwaves',
    description='Calculate wave-generated bottom orbital velocities from surface wave parameters',
    long_description=open('README.md').read(),
    setup_requires=['numpy>=1.7'],
    packages=find_packages(),
    package_data={'oceanwaves':['data/*']},
    entry_points={
        'console_scripts': [
            'spectralhs = oceanwaves.cmd.spectralhs:main',
            'ubspecclc3 = oceanwaves.cmd.ubspecclc3:main',
            'ubspecyv = oceanwaves.cmd.ubspecyv:main',
        ]
    },
)
