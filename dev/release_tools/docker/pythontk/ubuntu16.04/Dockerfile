FROM ubuntu:latest
MAINTAINER Nathan Dunfield <nathan@dunfield.info>

RUN apt-get update && apt-get install -y \
    openssh-server \
    x11-apps \
    python-tk \
    python-setuptools \
    python-pip \ 
    python3 \
    python3-setuptools \
    python3-pip \
    libglu1-mesa
RUN python -m pip install --upgrade pip setuptools virtualenv
RUN python3 -m pip install --upgrade pip setuptools virtualenv

# Setup sshd server
RUN mkdir /var/run/sshd
RUN echo 'root:snappy' | chpasswd
RUN echo '    PermitRootLogin yes' >> /etc/ssh/ssh_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

# Add non-root, but sudoable user.
RUN  adduser --quiet --shell /bin/bash --gecos "SnapPy user,101,," \ 
      --disabled-password snappy \
     && echo 'snappy:snappy' | chpasswd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

RUN mkdir /test
ADD test_pypi.py /test

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
