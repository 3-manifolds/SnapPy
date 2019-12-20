#!/bin/bash

set -e

cp -r /snappy .
cd snappy

echo
echo
echo "==================================================================="
echo "SageMath version"
sage --version

echo
echo
echo "==================================================================="
echo "Installing FXrays"
sage -pip install -U FXrays wheel

echo
echo
echo "==================================================================="
echo "Installing PLink"
sage -pip install -U https://github.com/3-manifolds/PLink/archive/master.zip

echo
echo
echo "==================================================================="
echo "Installing snappy_manifolds"
sage -pip install -U https://bitbucket.org/t3m/snappy_manifolds/get/tip.zip

echo
echo
echo "==================================================================="
echo "Installing spherogram"
sage -pip install -U https://bitbucket.org/t3m/spherogram/get/tip.zip

echo
echo
echo "==================================================================="
echo "Building SnapPy"
sage -python setup.py build

echo
echo
echo "==================================================================="
echo "Testing SnapPy"
sage -python setup.py test

echo
echo
echo "==================================================================="
echo "Building SnapPy docs"
sage -python setup.py build_docs
