=========================================
Validating our Linux install instructions
=========================================

This directory contains a suite of Docker specifications that mirror
(as of July 2019) our Linux install instructions to test that the
latter actualy, like, work.  Typical usage::

  docker build --tag=snappytest:ubuntu ubuntu
  
