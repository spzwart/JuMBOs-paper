#!/usr/bin/bash
## csh -f run_Pl.csh 0.25 -2.0 0.3 1
set rvir = $1
set x = $2
set mmin = $3
set nrun = $4

set alpha = 1.2
set Nstars = 2500
set NJupiters = 600
#set NJupiters = 5000

set t = 1
set dt = 0.1

set filename = "PlN"$Nstars"n"$NJupiters"x"$alpha"m"$mmin"ffc_R"$rvir"pcQ05_R"$nrun

mkdir $filename
cd $filename
/usr/bin/python3 ../../src/run_cluster.py --dt $dt -t $t -F 1.6 -R $rvir --Nstars $Nstars --NJuMBOs $NJupiters --mmin $mmin -x $x --model singletons > out.data &

