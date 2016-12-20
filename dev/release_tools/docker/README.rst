Repeatable Linux testing with Docker
====================================

Instructions are Docker 1.12 on an OS X host. Starting Docker via::

  open -a Docker

a docker "image" is a basically a powered-off VM state; a "container"
is a possibly running VM instantiated from an image::

  docker images     # Lists all images.
  docker ps -a      # All containers, running or not.

Images are deleted with "docker rmi", and containers with "docker rm".  

Images are built from "Dockerfiles" or pulled from an external
repository.  Key for use are the two "manylinux" containers
associated to::

  https://github.com/pypa/manylinux

To fetch them do::

  docker pull quay.io/pypa/manylinux1_x86_64
  docker pull quay.io/pypa/manylinux1_i686

To start one of these containers and open a shell on it do::

  docker run -i -t --name=test quay.io/pypa/manylinux1_x86_64 /bin/bash
  
For a second example, to build an image with Ubuntu 12.04 LTS, sshd,
python, tk and OpenGL do::

  docker build --tag pythontk:12.04   pythontk/ubuntu12.04

Now start it::

  docker run -d -P --name=test pythontk:12.04
  
and ssh in and test X11::

  docker ps     # See what port has been connected to port 22
  ssh -X root@localhost -p 32768  # password is "snappy"
  xclock &

To test the binary install of snappy, do::

  easy_install -U -f http://snappy.computop.org/get snappy
  python -m snappy.test

Now stop the container and destory it::

  docker stop test
  docker rm test

The end.


Images
======

We use two kinds of images; those for building binary eggs/wheels, and
those for testing.  The images for building are derived from the
manylinux images::

  docker build --tag many64 manylinux/64
  docker -i -t --name build64 many64





Testing binary builds::

  docker build --tag pythontk:12.04   pythontk/ubuntu12.04
  docker build --tag pythontk:14.04   pythontk/ubuntu14.04
  docker build --tag pythontk:centos6   pythontk/centos6

Current status:
============

Trying to get centos ssh to work, see

http://stackoverflow.com/questions/18173889/cannot-access-centos-sshd-on-docker


Rough version of build script
=============================


    cd build
    hg clone https://bitbucket.org/t3m/fxrays
    hg clone https://bitbucket.org/t3m/cypari
    hg clone https://bitbucket.org/t3m/plink
    hg clone https://bitbucket.org/t3m/spherogram
    hg clone https://bitbucket.org/t3m/snappy

    cd fxrays
    py27 setup.py sdist
    for py in py26, py27, py33, py34, py35
    py setup.py install
    py setup.py bdist_wheel
    py -m fxrays.test

    cd ../plink
    py27 setup.py sdist
    for py in py26, py27, py33, py34, py35:
    py setup.py install

    cd ../spherogram
    py27 setup.py sdist
    for py in py26, py27, py33, py34, py35:
    py setup.py install
    py setup.py bdist_wheel
    
    
    
    

    
