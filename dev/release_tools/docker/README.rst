Repeatable Linux testing with Docker
====================================

Instructions are for boot2docker VirtualBox VM on OS X host. Starting
the main VM::

  boot2docker start

Then issue all docker commands in a subshell where you have executed::
  
  bash
  $(boot2docker shellinit)

A docker "image" is a basically a powered-off VM state; a "container"
is a possibly running VM instantiated from an image::

  docker images     # Lists all images.
  docker ps -a      # All containers, running or not.

Images are deleted with "docker rmi", and containers with "docker rm".  

Images are built from "Dockerfiles". For example, to build an image
with Ubuntu 12.04 LTS, sshd, python, tk and OpenGL do::

  docker build --tag pythontk:12.04   pythontk/ubuntu12.04

Now start it::

  docker run -d -P --name=test pythontk:12.04
  
and ssh in and test X11::

  docker ps     # See what port has been connected to port 22
  boot2docker ip  # IP of VirtualBox VM
  ssh -X root@192.168.59.103 -p 49157  # password is "snappy"
  xclock &

To test the binary install of snappy, do::

  easy_install -U -f http://snappy.computop.org/get snappy
  python -m snappy.test

Now stop the container and destory it::

  docker stop test
  docker rm test

The end.  
