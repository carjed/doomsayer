from setuptools import setup, find_packages
setup(
    name="doomsayer",
    version="0.1",
    packages=find_packages(),
    install_requires=['cyvcf2', 'pyvcf', 'numpy', 'Biopython', 'scikit-learn'],

    # metadata for upload to PyPI
    author="Jedidiah Carlson",
    author_email="jed.e.carlson@gmail.com",
    description="Genome sequencing outlier detection using NMF",
    license="MIT",
    keywords="NMF genome sequencing",
    url="https://github.com/carjed/doomsayer",
)
