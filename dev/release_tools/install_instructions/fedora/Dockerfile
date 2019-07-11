FROM fedora:latest
LABEL maintainer="Nathan Dunfield <nathan@dunfield.info>"

RUN yum -y update; yum -y install sudo
RUN  adduser --shell /bin/bash snappy \
     && echo "snappy ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/01-snappy \
     && chmod 0440 /etc/sudoers.d/01-snappy

# Command from install instructions

USER snappy
RUN sudo yum -y install python3-tkinter python3-pip
RUN python3 -m pip install --upgrade --user snappy
RUN python3 -m snappy.test 
    
WORKDIR /home/snappy
CMD ["/bin/bash"]