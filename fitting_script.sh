#!/bin/bash

### fit paramaeters for a give simulation
#
#
make data-cleanup
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
      ./ovarian_follicle "./config/PhysiCell_settings.xml" "$k_oocyte" "$k_granulosa" "$k_basement"
    done
  done
done
#./ovarian_follicle "./config/PhysiCell_settings.xml" k_oocyte k_granulosa k_basement
