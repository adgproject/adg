# -*- coding: utf-8 -*-
"""Installation script for the Automated Diagaram Generator.

You can install the program either by running
    pip install <folder>
or
    python setup.py install
"""


import sys
from setuptools import setup
import adg

main_dependencies = [
    "setuptools"
]

for dep in main_dependencies:
    try:
        __import__(dep)
    except ImportError:
        print(
            "Error: You do not have %s installed, please\n"
            "       install it. For example doing\n"
            "\n"
            "       pip install %s\n" % (dep, dep)
        )
        sys.exit(1)


setup(
    name='adg',
    version=adg.__version__,
    maintainer='Pierre Arthuis',
    maintainer_email='pierre.arthuis@pm.me',
    author=adg.__author__,
    author_email=adg.__email__,
    license=adg.__license__,
    url='',
    install_requires=[
        "networkx>=2.0",
        "numpy",
    ],
    python_requires='>=2.7.1',
    extras_require=dict(
        # List additional groups of dependencies here (e.g. development
        # dependencies). You can install these using the following syntax:
        # $ pip install -e .[develop]
        develop=[
            'pytest',
            'pytest-cov',
            'sphinx',
            'sphinx_rtd_theme',
        ]
    ),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    description='A powerful diagram generator and evaluator for many-body '
                'formalisms in physics and chemistry',
    long_description='',
    keywords=[
        'physics',
        'chemistry',
        'theory',
        'many-body',
        'diagram',
        'cli',
        'batch'
    ],
    packages=[
        "adg",
    ],
    data_files=[

        ("share/doc/adg/", [
            "README.md",
        ]),

        ("share/man/man1", [
            "doc/build/man/adg.1",
        ]),

    ],
    entry_points=dict(
        console_scripts=[
            'adg=adg.main:main'
        ]
    ),
    platforms=['linux', 'osx'],
)
