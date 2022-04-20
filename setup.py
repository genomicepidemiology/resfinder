import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="resfinder",
    version="2.4.0",
    author="Center for Genomic Epidemiology",
    author_email="food-cgehelp@dtu.dk",
    description=("ResFinder identifies acquired genes and/or finds chromosomal "
                 "mutations mediating antimicrobial resistance in total or "
                 "partial DNA sequence of bacteria."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/genomicepidemiology/resfinder",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 1 - Planning",
    ],
    python_requires='>=3.10',
    keywords='bioinformatics antimicrobial resistance',
    project_urls={
        'Repository': 'https://bitbucket.org/genomicepidemiology/resfinder',
        'Bug Tracker': ('https://bitbucket.org/genomicepidemiology/resfinder'
                        '/issues'),
        'CGE Website': 'http://genomicepidemiology.org',
        'CGE Webtools': 'http://genomicepidemiology.org/services',
    },

)
