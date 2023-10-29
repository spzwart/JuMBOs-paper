set t=$1
set dir = $2
#set outfile = ../data_reduced/$dir.$t.csv

#echo "type,key,k1,k2,M,m,a,e,i,x,y,z,vx,vy,vz" > $outfile

set i = 1
while ( $i <= 10 )
    set dir1 = $dir$i/jumbos_i0000.amuse
    set dir2 = $dir$i/jumbos_i$t.amuse
    python ../data_reduced/find_collision_products.py --f1 $dir1 --f2 $dir2 | grep "Number"
   @ i++
end


