#!/usr/bin/env bash

for dir in ./*.ms
do
    dir=${dir%*/}
    echo ${dir##*/}
    tar -zcvf ${dir##*/}.tar.gz ${dir##*/}
done

