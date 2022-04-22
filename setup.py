import os.path
import setuptools


VERSION_FILE = "VERSION"
README_FILE = "README.md"
CHANGELOG_FILE = "CHANGELOG.md"


def read_property(filename, dirs=""):
    filepath = os.path.join(os.path.dirname(__file__), dirs, filename)
    with open(filepath, "r", encoding="utf-8") as fh:
        return fh.read()


long_description = ("{readme}\n----------\n{change}"
                    .format(readme=read_property(README_FILE),
                            change=read_property(CHANGELOG_FILE)))

setuptools.setup(
    name="resfinder",
    version=read_property(VERSION_FILE, dirs="src/resfinder/").strip(),
    author="Center for Genomic Epidemiology",
    author_email="food-cgehelp@dtu.dk",
    description=("ResFinder identifies acquired genes and/or finds chromosomal "
                 "mutations mediating antimicrobial resistance in total or "
                 "partial DNA sequence of bacteria."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/genomicepidemiology/resfinder",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 1 - Planning",
    ],
    packages=setuptools.find_packages(where='src'),
    package_dir={"": "src"},
    python_requires='>=3.10',
    install_requires=[
        'cgelib',
        'cgecore==1.5.6',
        'tabulate',
        'biopython',
        'pandas'
    ],
    include_package_data=True,
    keywords='bioinformatics antimicrobial resistance',
    project_urls={
        'Repository': 'https://bitbucket.org/genomicepidemiology/resfinder',
        'Bug Tracker': ('https://bitbucket.org/genomicepidemiology/resfinder'
                        '/issues'),
        'CGE Website': 'http://genomicepidemiology.org',
        'CGE Webtools': 'http://genomicepidemiology.org/services',
    },

)
