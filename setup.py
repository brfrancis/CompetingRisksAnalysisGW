import os

from setuptools import setup

setup(
    name="crag",
    author="Ben Francis",
    author_email="b.r.francis@liverpool.ac.uk",
    description="Competing Risk Analysis Genome-wide",
    keywords="competing risks survival analysis statistics data analysis",
    url="https://github.com/brfrancis/CompetingRisksAnalysisGW",
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
