FROM python:3.7-stretch
MAINTAINER Nathan Dunfield <nathan@dunfield.info>

RUN apt-get update && apt-get install -y \
    python3-tk \
    libglu1-mesa-dev

RUN python3 -m pip install -U pip setuptools wheel ipython
RUN python3 -m pip install cypari snappy_manifolds
RUN python3 -m pip install --no-binary :all: snappy
RUN python3 -m snappy.test

WORKDIR /
CMD ["/bin/bash"]
