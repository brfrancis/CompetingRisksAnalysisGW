import os

from setuptools import setup

setup(
    name="crag",
    author="Ben Francis",
    author_email="b.r.francis@liverpool.ac.uk",
    description="Competing Risk Analysis Genome-wide",
    keywords="competing risks survival analysis statistics data analysis",
    url="https://github.com/brfrancis/CompetingRisksAnalysisGW",
    packages=['crag'],
    long_description=readme_rst,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering",
        ],
    install_requires=[
      "subprocess",
      "sys",
      "math",
      "csv",
      "pandas",
      "gzip",
      "time",
      "argparse",
      "warnings",
        "numpy",
        "lifelines",
        "pandas>=0.18",
    ],
)
subprocess, sys, math, lifelines, csv, numpy, pandas, gzip, time, argparse and warnings
