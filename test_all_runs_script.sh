#!/bin/bash
# make clean
rm -rf ./output
mkdir ./output
touch ./output/empty.txt
rm -rf ./all_runs
mkdir ./all_runs
touch ./all_runs/empty.txt
# make
./ovarian_follicle "./config/EG-withCells_0p1.xml" >>./output/output_log
mv  ./output/ ./all_runs/EG-ouput
mkdir ./output
./ovarian_follicle "./config/EGandGLY-withCells_0p01.xml" >> ./output/output_log
mv  ./output/ ./all_runs/EGandGLY-ouput
mkdir ./output
./ovarian_follicle "./config/GLY-noCells_0p1.xml" >> ./output/output_log 
mv  ./output/ ./all_runs/GLY-ouput
mkdir ./output
./ovarian_follicle "./config/HM05-withCells_0p1.xml" >> ./output/output_log 
mv  ./output/ ./all_runs/HM05-ouput
mkdir ./output
./ovarian_follicle "./config/HM1-withCells_0p1.xml" >> ./output/output_log 
mv  ./output/ ./all_runs/HM1-ouput
mkdir ./output
./ovarian_follicle "./config/HM2-withCells_0p1.xml" >> ./output/output_log 
mv  ./output/ ./all_runs/HM2-ouput
mkdir ./output
./ovarian_follicle "./config/HM5-withCells_0p1.xml" >> ./output/output_log 
mv  ./output/ ./all_runs/HM5-ouput
mkdir ./output
