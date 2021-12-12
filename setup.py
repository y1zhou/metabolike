from setuptools import find_packages, setup

from metabolike import __version__

setup(
    name="metabolike",
    verion=__version__,
    author="Yi Zhou",
    url="https://github.com/y1zhou/metabolike",
    packages=find_packages(),
    keywords=["bioinformatics", "metabolic-pathway", "graph"],
    python_requires=">=3.6",
    install_requires=["lark>=0.11.3,<1.0.0", "pandas", "numpy", "pyspark"],
)
