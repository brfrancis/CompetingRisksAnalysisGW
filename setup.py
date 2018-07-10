import os

from setuptools import setup, find_packages

datloc="data/*"
bdatloc=str.encode(datloc)



def filepath(fname):
    return os.path.join(os.path.dirname(__file__), fname)

readme_md = filepath('README.md')

try:
    import pypandoc
    readme_rst = pypandoc.convert_file(readme_md, 'rst')
except(ImportError):
    readme_rst = open(readme_md,encoding='utf-8').read()

setup(
    name="crag",
    version="1.0.5",
    author="Ben Francis",
    author_email="b.r.francis@liverpool.ac.uk",
    description="Competing Risk Analysis Genome-wide",
    license="GPLv3",
    long_description=readme_rst,
    long_description_content_type="text/markdown",
    keywords="competing risks survival analysis statistics data analysis",
    url="https://github.com/brfrancis/CompetingRisksAnalysisGW",
    packages=['crag'],
    classifiers=(
      "Programming Language :: Python :: 2",
      "Programming Language :: Python :: 3",
      "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
      "Operating System :: OS Independent",
    ),
    install_requires=[
      "pandas",
      "argparse",
      "numpy",
      "lifelines"
    ],
    package_data={
      'crag': [
        "data/*",
      ]
    },
    scripts=['bin/crag_gwas.py'],
)
