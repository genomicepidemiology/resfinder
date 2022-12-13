Development documentation
=============

This document is intended for developers working on the ResFinder application. It consists of suggestions on how to set up a development environment and guidelines for doing releases of ResFinder.

## Setup development environment

### 1. Install PDM - Python Development Master

In order to test and build ResFinder as described here, you will need PDM.

*"PDM is a modern Python package manager with PEP 582 support. It installs and manages packages in a similar way to npm that doesn't need to create a virtualenv at all!"* - **PDM website**

* Install PDM as described on the [PDM website](https://daobook.github.io/pdm/).

### 2. Install BLAST and KMA

### 3. Clone ResFinder repository

### 4. Clone Databases and index them

For the tests to work, you need to set the environment variables for the database locations.

```bash

export CGE_RESFINDER_RESGENE_DB="/path/to/resfinder_db/"
export CGE_RESFINDER_RESPOINT_DB="/path/to/pointfinder_db/"

```

### 5. Setup PDM with ResFinder dependencies

```bash

# Go to ResFinder root directory
cd /path/to/resfinder/

# Install all python depencies inside the directory (resfinder/__pypackages__)
pdm install

```

## Running and testing ResFinder

```bash

# Go to ResFinder root directory
cd /path/to/resfinder/

# Run ResFinder
pdm run resfinder -h

# Run ResFinder tests
pdm run tests

```

## Creating / Edit tests

All tests are written as doctests. The doctests are stored in markdown formatted files in the `tests` directory. The `tests` directory mirrors the structure of the `src` directory. Tests are written into files named after a corresponding file in the `src` directory tree, but prefixed with `test_`. A test is written into the test file that corresponds to the file in which the tested code resides.

*Example*

**Code**: `src/resfinder/cge/output/gene_result.py`

**Test file**: `tests/resfinder/cge/output/test_gene_result.md`

If a `*.py` file doesn't have a corresponding test file you need to create it if you need to write tests for the code. Remember to prefix the name with `test_`, as this will automatically include the tests when running `pdm run tests`. Use existing test files to get an idea of how these should be formatted.

**Note**: The configuration for recognizing tests are set in the file `pytest.ini` in the root directory.


## Deploy

1. Change version number in `src/resfinder/__init__.py`.
2. Make sure CHANGELOG.md is up to date and add the current date of release.
3. Build package:

```bash

pdm build

```

4. Upload package to pypi:

```bash

twine upload dist/*

```

### Deploy docker image

**ToDo**