FROM quay.io/pypa/manylinux1_x86_64:latest
MAINTAINER Nathan Dunfield <nathan@dunfield.info>

RUN yum install -y nano tk mesa-libGLU-devel
RUN ln -s /opt/python/cp26-cp26mu/bin/python /bin/py26
RUN ln -s /opt/python/cp27-cp27mu/bin/python /bin/py27
RUN ln -s /opt/python/cp33-cp33m/bin/python  /bin/py33
RUN ln -s /opt/python/cp34-cp34m/bin/python  /bin/py34
RUN ln -s /opt/python/cp35-cp35m/bin/python  /bin/py35
RUN py27 -m pip install mercurial
RUN ln -s /opt/python/cp27-cp27mu/bin/hg /bin/hg
RUN mkdir /build
RUN /opt/python/cp26-cp26mu/bin/pip install cython
RUN py27 -m pip install cython
RUN py33 -m pip install cython
RUN py34 -m pip install cython
RUN py35 -m pip install cython

WORKDIR /build
CMD ["/bin/bash"]