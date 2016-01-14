FROM centos:6

RUN echo 'root:snappy' | chpasswd
RUN yum update -y && yum install -y \
    openssh-server \
    xorg-x11-apps \
    xorg-x11-xauth \
    tkinter \
    python-setuptools \
    mesa-libGLU

ADD ssh_host_* /etc/ssh/

EXPOSE 22
CMD    ["/usr/sbin/sshd", "-D"]
