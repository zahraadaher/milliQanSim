#!/bin/bash
#create local build of sim, configured to generate cosmic muons
. buildsetup.sh
cp CMakeLists.txt.default CMakeLists.txt
cp include/default/* include/
cp include/defaultSlab/* include/
cp src/default/* src/
cp src/defaultSlab/* src/

# Use versioned shielding file for Geant4 version 11
# Copies all files with .Geant4.* suffix to the same name without the suffix

for f in src/*.Geant4.11*; do
    base="${f%%.Geant4.11*}"   # remove .Geant4.* from filename
    echo "Copying $f -> ${base}"
    cp "$f" "${base}"
done

cd build
cp -r ../inputData/config .
cmake ../
make -j8
./MilliQan ../runMac/specmuonFullPropagated_no_vis.mac
