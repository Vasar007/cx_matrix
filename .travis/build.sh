#!/bin/bash

set -x

mkdir build
cd build

cmake .. -DCMAKE_BUILD_TYPE=$BUILD_CONFIGURATION -DTARGET_CPU=$TARGET_CPU || exit 1
cmake --build . || exit 1
ctest --output-on-failure || exit 1
