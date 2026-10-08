#!/bin/bash
#create local build of sim, configured to generate signal mcps for the slab detector
cp inputData/config/particlesMCP.ini inputData/config/particles.ini
. buildsetup.sh
cp MilliQan.cc.BeamGen MilliQan.cc
cp CMakeLists.txt.default CMakeLists.txt
cp include/default/* include/
cp src/default/*.cc src/
cp include/defaultSlab/* include/
cp src/defaultSlab/* src/

cd build
cmake ../
make -j8
