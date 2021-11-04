FROM quay.io/pypa/manylinux2010_x86_64:latest
MAINTAINER Nathan Dunfield <nathan@dunfield.info>

RUN yum install -y nano tk mesa-libGLU-devel openssh-clients wget
RUN ln -s /opt/python/cp36-cp36m/bin/python  /bin/py36
RUN ln -s /opt/python/cp37-cp37m/bin/python  /bin/py37
RUN ln -s /opt/python/cp38-cp38/bin/python   /bin/py38
RUN ln -s /opt/python/cp39-cp39/bin/python   /bin/py39
RUN ln -s /opt/python/cp310-cp310/bin/python /bin/py310
RUN py38 -m pip install --no-warn-script-location twine

# Downgrade auditwheel by one notch so that it doesn't see the extra
# libs in the Togl binaries. The cause is this::
#
# https://github.com/pypa/auditwheel/pull/95
#
# and the problem is that "repair" fails. 
# RUN py36 -m pip install --no-warn-script-location auditwheel==1.8.0

RUN mkdir /build
RUN py36 -m pip install --upgrade --no-warn-script-location \
    pip setuptools cython sphinx decorator future ipython networkx
RUN py37 -m pip install --upgrade --no-warn-script-location \
    pip setuptools cython sphinx decorator future ipython networkx
RUN py38 -m pip install --upgrade --no-warn-script-location \
    pip setuptools cython sphinx decorator future ipython networkx
RUN py39 -m pip install --upgrade --no-warn-script-location \
    pip setuptools cython sphinx decorator future ipython networkx
RUN py310 -m pip install --upgrade --no-warn-script-location \
    pip setuptools cython sphinx decorator future ipython networkx


ADD bin /build/bin
RUN chmod +x /build/bin/*
WORKDIR /build
CMD ["/bin/bash"]