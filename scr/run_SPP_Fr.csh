#!/usr/bin/bash
set rvir = $1
set amin = $2
set amax = $3
set nrun = $4

set filename = "FrN2500n300spp_a"$amin+$amax"R"$rvir"pcQ05_R"$nrun

mkdir $filename
cd $filename
/usr/bin/python3 ../../src/run_cluster.py --a1 $amin --a2 $amax -F 1.6 -R $rvir --NJuMBOs 300 --model circum_stellar > out.data &

