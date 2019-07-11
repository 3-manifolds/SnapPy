FROM archlinux/base
LABEL maintainer="Nathan Dunfield <nathan@dunfield.info>"

RUN  pacman --sync --refresh
RUN  pacman --sync --noconfirm sudo
RUN  useradd --create-home snappy \
     && echo "snappy ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-snappy \
     && chmod 0440 /etc/sudoers.d/01-snappy

# Commands from install instructions.  Note that in arch, python is
# Python 3, not Python 2.

USER snappy
RUN sudo pacman -Sy --noconfirm python-pip tk
RUN python -m pip install --upgrade --user snappy
RUN python -m snappy.test 
    
WORKDIR /home/snappy
CMD ["/bin/bash"]