# coding: utf-8
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='interface_stability',
    version='0.1',
    description='Python project for analyze interface stability',
    long_description=long_description,
    url='https://github.com/mogroupumd/interface_stability',
    author='Yifei Mo Research Group at UMD',
    author_email='yizhou.zhu@gmail.com',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6'
    ],
    keywords=[],
    packages=find_packages(),
    install_requires=['pymatgen>=2018.5.3','argparse'],
    extras_require={},
    package_data={},
    data_files=[],
    entry_points={
        'console_scripts': [
            'phase_stability=interface_stability.scripts.phase_stability:main',
            'pseudo_binary=interface_stability.scripts.pseudo_binary:main',
        ],
    },
    project_urls={
        'Bug Reports': 'https://github.com/mogroupumd/interface_stability/issues',
        'Source': 'https://github.com/mogroupumd/interface_stability',
    },
)