# Development documentation

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

# Some tests requires the following environmental variables to be set
export CGE_BLASTN="/path/to/blastn"
export CGE_RESFINDER_RESGENE_DB="/path/to/resfinder_db"
export CGE_RESFINDER_RESPOINT_DB="/path/to/pointfinder_db"
export CGE_DISINFINDER_DB="/path/to/disinfinder_db"

# Run ResFinder tests
pdm run test

```

## Creating / Edit tests

All tests are written as doctests. The doctests are stored in markdown formatted files in the `tests` directory. The `tests` directory mirrors the structure of the `src` directory. Tests are written into files named after a corresponding file in the `src` directory tree, but prefixed with `test_`. A test is written into the test file that corresponds to the file in which the tested code resides.

*Example*:

**Code**: `src/resfinder/cge/output/gene_result.py`

**Test file**: `tests/resfinder/cge/output/test_gene_result.md`

If a `*.py` file doesn't have a corresponding test file you need to create it if you need to write tests for the code. Remember to prefix the name with `test_`, as this will automatically include the tests when running `pdm run tests`. Use existing test files to get an idea of how these should be formatted.

**Note**: The configuration for recognizing tests are set in the file `pytest.ini` in the root directory.

## Deploy

0. If you plan to also make a Docker image release go to section [Deploy docker image](#Deploy docker image) and do steps 1-4.
1. Change version number in `src/resfinder/__init__.py`.
2. Make sure CHANGELOG.md is up to date and add the current date of release.
3. Push changes to repository.
4. Tag the commit you just pushed with the correct version number.
5. Build package:

    ```bash

    pdm build

    ```

6. Upload package to pypi:

    ```bash

    twine upload dist/*

    ```

### Docker image

**Note**: We do not guarentee that all ResFinder releases has a coresponding Docker image release. However, we should strive toward making as many as possible. At least Docker image releases should be done whenever significant changes has been released.

A docker image contains both the ResFinder software and the three databases:
[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/)
[PointFinder](https://bitbucket.org/genomicepidemiology/pointfinder_db/)
[DisinFinder](https://bitbucket.org/genomicepidemiology/disinfinder_db/)

#### Versioning

Versioning is done so that a specific version of a Docker image can always track exactely which version of the ResFinder sofware is included and exactely which databases.

* A Docker release version should match the ResFinder version.
* All databases included in a Docker image release should be tagged with a version number for that database release.
* The commits that correspond to the database versions used should be tagged with `resfinder-VERSION`, where `VERSION` matches the ResFinder and therfore also the Docker image version. Hence, all database versions/commits included will contain at least two tags: a database version number and `resfinder-VERSION`.

#### Deploy docker image

1. Make sure each database you want to include is tagged with a version number.
2. Make sure you have a released (versioned) ResFinder repo.
3. Tag each database with `resfinder-VERSION`, where `VERSION` matches the ResFinder version number (ex.: `resfinder-4.2.3`).
4. Update the Dockerfile:
    * (optional) Bump the KMA version (no major version change).
    * Change RESFINDER_VERSION environment variable.
5. Build Docker image. Note `<VERSION>` should be replaced with ResFinder version number:

    ```bash

    # Go to ResFinder root directory
    cd /path/to/resfinder/

    # Build image - No cache is used as you risc building with an outdated cached database otherwise.
    docker build --no-cache -t genomicepidemiology/resfinder:<VERSION> .

    ```

6. Push Docker image:

    ```bash

    docker push genomicepidemiology/resfinder:<VERSION>

    ```

7. Update latest tag:

    ```bash

    # Get image id, replace with <ID> later
    $ docker images
    REPOSITORY                     TAG        IMAGE ID       CREATED        SIZE
    cgetools_front-web             latest     cf74e74364c9   2 months ago   392MB
    cgetools_front-celery_worker   latest     6524fa33a75a   2 months ago   392MB
    redis                          7-alpine   39267c75a230   2 months ago   28.1MB
    rf_test                        latest     6692e77b787f   3 months ago   760MB
    
    $ docker tag <ID> genomicepidemiology/resfinder:latest
    $ docker push genomicepidemiology/resfinder:latest

    ```
