FROM python:3.7-rc-stretch
MAINTAINER Nathan Dunfield <nathan@dunfield.info>

RUN apt-get update && apt-get install -y \
    python3-tk \
    libglu1-mesa-dev

RUN python3 -m pip install --upgrade pip setuptools wheel cython sphinx decorator future ipython networkx

RUN hg clone https://bitbucket.org/t3m/cypari /opt/cypari
RUN hg clone https://bitbucket.org/t3m/snappy /opt/snappy

WORKDIR /opt/cypari/Version2
RUN python3 setup.py build install

RUN python3 -m pip install hg+https://bitbucket.org/t3m/plink
RUN python3 -m pip install hg+https://bitbucket.org/t3m/spherogram

WORKDIR /opt/snappy
RUN python3 setup.py build -j 4
RUN python3 setup.py pip_install
RUN python3 -m snappy.test

WORKDIR /
CMD ["/bin/bash"]
