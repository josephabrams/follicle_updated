#!/bin/bash
# note in main.cpp run time is specified staticlly to overwrite xml
# make clean
rm -rf ./output
mkdir ./output
touch ./output/empty.txt
mkdir ./fitting
# make
### fit paramaeters for a give simulation
#
#
k_oocyte="0"
k_granulosa="0"
k_basement="0"
for i in {0..9}
do
  k_oocyte="0.$i"
  for j in {0..9}
  do
    k_granulosa="0.$j"
    for k in {0..9} 
    do
      k_basement="0.$k"
      #echo "$k_oocyte, $k_granulosa, and $k_basement" 
      ./ovarian_follicle "./config/EG-withCells_0p1.xml" "$k_oocyte" "$k_granulosa" "$k_basement"
      
    done
  done
done

mv  ./output/ ./fitting/EG-fitting
mkdir ./output
#./ovarian_follicle "./config/PhysiCell_settings.xml" k_oocyte k_granulosa k_basement
