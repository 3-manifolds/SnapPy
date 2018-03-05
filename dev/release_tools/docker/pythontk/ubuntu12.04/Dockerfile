FROM ubuntu:12.04
MAINTAINER Nathan Dunfield <nathan@dunfield.info>

RUN apt-get update && apt-get install -y \
    openssh-server \
    x11-apps \
    python-tk \
    python-setuptools \
    libglu1-mesa
RUN mkdir /var/run/sshd
RUN echo 'root:snappy' | chpasswd

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
