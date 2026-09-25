#!/bin/sh

for i in 2.4603
do
    mkdir acc$i
    cd acc$i
    cp ../acc2.46/GenMoire.f90 .
    echo $i
    sed -i -e "s/2.46/$i/g" GenMoire.f90 > GenMoireNew.f90
    gfortran GenMoireNew.f90
    ./a.out
    cd ..
done
