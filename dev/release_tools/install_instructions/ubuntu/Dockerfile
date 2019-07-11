FROM ubuntu:18.04
LABEL maintainer="Nathan Dunfield <nathan@dunfield.info>"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y sudo tzdata

# Add a nonroot user.

RUN  adduser --quiet --shell /bin/bash --gecos "SnapPy user,101,," \
               --disabled-password snappy \
     && echo "snappy ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-snappy \
     && chmod 0440 /etc/sudoers.d/01-snappy
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Command from install instructions

USER snappy
RUN sudo apt-get install -y python3-tk python3-pip
RUN python3 -m pip install --upgrade --user snappy
RUN python3 -m snappy.test


WORKDIR /home/snappy
CMD ["/bin/bash"]