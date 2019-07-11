FROM ubuntu:18.04
LABEL maintainer="Nathan Dunfield <nathan@dunfield.info>"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y sudo tzdata wget

# Add a nonroot user.

RUN  adduser --quiet --shell /bin/bash --gecos "SnapPy user,101,," \
               --disabled-password snappy \
     && echo "snappy ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-snappy \
     && chmod 0440 /etc/sudoers.d/01-snappy
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Command from install instructions

USER snappy
RUN sudo apt-get install -y python3-tk python3-dev python3-pip libglu1-mesa-dev
RUN python3 -m pip install --user plink
RUN python3 -m pip install --user cython

WORKDIR /home/snappy
RUN wget https://pypi.python.org/packages/source/s/snappy/snappy-2.6.tar.gz
RUN tar xfz snappy-2.6.tar.gz
WORKDIR /home/snappy/snappy-2.6
RUN python3 setup.py build
RUN python3 -m pip install --user .
RUN python3 -m snappy.test

CMD ["/bin/bash"]