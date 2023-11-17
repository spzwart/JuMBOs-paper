#!/usr/bin/bash
set rvir = $1
set amin = $2
set amax = $3
set nrun = $4

set filename = "PlN2500n300isf_a"$amin+$amax"R"$rvir"pcQ05_R"$nrun

mkdir $filename
cd $filename
/usr/bin/python3 ../../src/run_cluster.py --a1 $amin --a2 $amax -R $rvir --NJuMBOs 300 --model freefloaters > out.data &

