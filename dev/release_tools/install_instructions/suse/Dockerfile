FROM opensuse/leap:latest
LABEL maintainer="Nathan Dunfield <nathan@dunfield.info>"


RUN  zypper install -y sudo
RUN  useradd --create-home snappy \
     && echo "snappy ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-snappy \
     && chmod 0440 /etc/sudoers.d/01-snappy

# Command from install instructions

USER snappy
RUN sudo zypper install -y python3-tk python3-pip
RUN python3 -m pip install --upgrade --user snappy
RUN python3 -m snappy.test
   
WORKDIR /home/snappy
CMD ["/bin/bash"]