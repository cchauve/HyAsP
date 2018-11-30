#!/usr/bin/env python

from setuptools import setup

setup(
    name = 'HyAsP',
    version = '1.0.0',
    description = 'HyAsP (Hybrid Assember for Plasmids) determines plasmids in assemblies based on characteristics of contigs and the occurrences of plasmid genes.',
    author = 'Robert Mueller <romueller@techfak.uni-bielefeld.de>, Cedric Chauve <cedric_chauve@sfu.ca>',
    author_email = 'cedric_chauve@sfu.ca',
    url = 'https://github.com/cchauve/HyAsP',
    packages = ['HyAsP'],
    license = 'MIT',

    install_requires = [
        'biopython>=1.71',
        'numpy>=1.14.3',
        'pandas>=0.22.0'
    ],

    entry_points = {
        'console_scripts': [
            'hyasp=HyAsP.hyasp:main',
            'hyasp_pipeline=HyAsP.fastq_to_plasmids:main'
        ]
    }
)
