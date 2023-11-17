set t=$1
set dir = $2
set run = $3
set outfile = ../data_reduced/$dir.R$run.$t.csv

echo "type,key,k1,k2,M,m,a,e,i,x,y,z,vx,vy,vz" > $outfile

set i = 1
while ( $i <= $run )
    set dirname = $dir$i/jumbos_i$t.amuse
    /usr/bin/python3 ../src/find_multiples.py -f $dirname >>  $outfile
   @ i++
end


#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R2/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R3/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R4/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R5/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R6/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R7/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R8/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R9/jumbos_i$t.amuse >>  $outfile
#python ../src/find_multiples.py -f PlN2500n600cs_R1.0pcQ05_R10/jumbos_i$t.amuse >> $outfile

#cat PlN2500n600cs_R1.0pcQ05fm_t$t.out  | grep 'binary\|triple\|triplesing\|binbin' | grep -v 'N ' >> ../data_reduced/PlN2500n600cs_R1.0pcQ05fm_t$t.csv


# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R11/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R12/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R13/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R14/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R15/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R16/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R17/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R18/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R19/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R20/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out

# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R21/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R22/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R23/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R24/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R25/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R26/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R27/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R28/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R29/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R30/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out

# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R31/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R32/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R33/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R34/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R35/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R36/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R37/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R38/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R39/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
# python ../src/find_multiples.py -f PlN2500n600cs_R0.25pcQ05_R40/jumbos_i$t.amuse >>  PlN2500n600cs_R0.25pcQ05fm_t$t.out
