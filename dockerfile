FROM python:3.10-slim-bullseye

ENV INSIDE_DOCKER true
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -qq; \
    apt-get install --no-install-recommends -y -qq git \
    build-essential \
    wget \
    ncbi-blast+ \
    libz-dev \
    procps \
    ; \
    rm -rf /var/cache/apt/* /var/lib/apt/lists/*;

# Install KMA
RUN cd /usr/src/; \
    git clone --depth 1 -b 1.4.11 https://bitbucket.org/genomicepidemiology/kma.git; \
    cd kma && make; \
    mv kma /usr/bin/; \
    cd ..; \
    rm -rf kma/;

ENV RESFINDER_VERSION 4.3.3

# Install ResFinder
RUN pip install --no-cache-dir resfinder==$RESFINDER_VERSION

# Install databases
RUN cd /; \
    mkdir databases; \
    cd /databases/; \
    git clone -b resfinder-$RESFINDER_VERSION --depth 1 https://git@bitbucket.org/genomicepidemiology/resfinder_db.git; \
    git clone -b resfinder-$RESFINDER_VERSION --depth 1 https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git; \
    git clone -b resfinder-$RESFINDER_VERSION --depth 1 https://git@bitbucket.org/genomicepidemiology/disinfinder_db.git; \
    cd /databases/resfinder_db; \
    python INSTALL.py; \
    rm -rf .git; \
    cd /databases/pointfinder_db; \
    python INSTALL.py; \
    rm -rf .git; \
    cd /databases/disinfinder_db; \
    python INSTALL.py; \
    rm -rf .git;

# Setup environment
RUN cd /; \
    mkdir app;
WORKDIR /app

# Environmental variables
ENV CGE_RESFINDER_RESGENE_DB /databases/resfinder_db/
ENV CGE_RESFINDER_RESPOINT_DB /databases/pointfinder_db/
ENV CGE_DISINFINDER_DB /databases/disinfinder_db/

# Execute program when running the container
ENTRYPOINT ["python", "-m", "resfinder"]
