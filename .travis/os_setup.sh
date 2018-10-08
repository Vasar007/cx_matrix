#!/bin/bash

set -x

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt-get update -qq

    CMAKE_VERSION=3.8.2
    CMAKE_VERSION_DIR=v3.8

    CMAKE_OS=Linux-x86_64
    CMAKE_TAR=cmake-$CMAKE_VERSION-$CMAKE_OS.tar.gz
    CMAKE_URL=http://www.cmake.org/files/$CMAKE_VERSION_DIR/$CMAKE_TAR
    CMAKE_DIR=$(pwd)/cmake-$CMAKE_VERSION

    wget --quiet $CMAKE_URL
    mkdir -p $CMAKE_DIR
    tar --strip-components=1 -xzf $CMAKE_TAR -C $CMAKE_DIR
    export PATH=$CMAKE_DIR/bin:$PATH
fi