=================================================
Repeatable Linux building and testing with Docker
=================================================

Introduction
============

The instructions here are for Docker 1.12 on an OS X host assuming all
commands are executed from the directory containing this file.  First,
start Docker via::

  open -a Docker

A docker "image" is a basically a powered-off VM state; a "container"
is a possibly running VM instantiated from an image::

  docker images     # Lists all images.
  docker ps -a      # All containers, running or not.

Images are deleted with "docker rmi", and containers with "docker rm".  

Images are built from "Dockerfiles" or pulled from an external
repository.  Key for us are the "manylinux1" images associated to::

  https://github.com/pypa/manylinux

To fetch the 64-bit one do::

  docker pull quay.io/pypa/manylinux1_x86_64
  
To start this container and open a shell on it do::

  docker run -i -t --name=test quay.io/pypa/manylinux1_x86_64 /bin/bash

After logging out, the container automatically stops, but you probably
want to also delete it with::
  
  docker rm test


Building binaries: Manylinux image
==================================

Our Linux binaries are built using an image derived from the above
manylinux1 image via::

  docker build --tag=many64 manylinux
  
To start a container with this image and open a shell on it do::

  docker run -i -t --name=build64 many64

When you log out, while the container is stopped any changes you made
are saved.  You can restart the container via::

  docker start -i build64

and create addtional shells on a running container via::

  docker exec -i -t build64 /bin/bash

Once logged in, you will be in the "/build" directory, and there are
aliases "py27", "py34", etc. for the various versions of Python.
Additionally, "/build/bin" contains three scripts that make it easy to
build Linux binaries. Typical usage::

  bin/clone                # Gets latest sources from Bitbucket.
  bin/build FXrays py27    # Builds Python 2.7 eggs/wheels for FXrays.
  py27 -m FXrays.test -v   # Run doctests.
  bin/build FXrays         # Builds all supported eggs/wheels for FXrays.
  bin/collect              # Copies eggs/wheels into /build/dist.
  scp -r dist dunfield@thurston.math.illinois.edu:Dropbox/snappy-release

Note that the "collect" script runs "auditwheels" to ensure that
our wheels are "manylinux1" compliant.


Images for testing
==================

**Major Problem**: As of 2016/12, XQuartz on OS X is too broken for
OpenGL graphics to work over X11, rendering this section moot.  So
while e.g. xclock works fine, the full glory of SnapPy is
thwarted. For this reason, Linux testing is currently done in a
VMWare Fusion VM.

**Different Approach**: Tunnelling X11 over ssh is excessive, there is
easier way, though it doesn't help with the previous
problem. Following::

https://fredrikaverpil.github.io/2016/07/31/docker-for-mac-and-gui-applications/

one sets XQuartz to allow incoming connections, and then on the macOS
side do::

  xhost +130.126.111.217    # = thurston.math.illinois.edu
  
Then open a shell on a Docker container and do::

  export DISPLAY=130.126.111.217:0
  xclock &

To build an image with Ubuntu 14.04 LTS, sshd, python 2 and 3, tk and
OpenGL do::

  docker build --tag pythontk:14.04   pythontk/ubuntu14.04

Now start it::

  docker run -d -P --name=test14 pythontk:14.04
  
and ssh in and test X11::

  docker ps     # See what port has been connected to port 22
  ssh -X root@localhost -p 32768  # password is "snappy"
  xclock &

There is a script for automatically downloading and testing one of our
packages, which you can use like this::

  cd /test
  python test_pypi.py -p FXrays

For more options, do::

  python test_pypi.py -h

To test with Python 3 rather than Python 2, just run "test_pypi.py"
with "python3".  To stop the container and destory it, do::

  docker stop test
  docker rm test

The end.
 
