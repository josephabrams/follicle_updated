#!/bin/bash

./ovarian_follicle "./config/EG-withCells_0p0001.xml" >> ./output/output_log 
mv ./output/ ./EG-withCells_0p0001
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p001.xml" >> ./output/output_log 
mv ./output/ ./EG-withCells_0p001
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p01.xml" >> ./output/output_log 
mv ./output/ ./EG-withCells_0p01
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p02.xml" >> ./output/output_log
mv ./output/ ./EG-withCells_0p02
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p05.xml" >> ./output/output_log
mv ./output/ ./EG-withCells_0p05
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p1.xml" >> ./output/output_log
mv ./output/ ./EG-withCells_0p1
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p2.xml" >> ./output/output_log
mv ./output/ ./EG-withCells_0p2
mkdir ./output

./ovarian_follicle "./config/EG-withCells_0p5.xml" >> ./output/output_log
mv ./output/ ./EG-withCells_0p5
mkdir ./output

# ./ovarian_follicle "./config/EGandGLY-withCells_0p0001.xml" >> ./output/output_log 
# mv ./output/ ./EGandGLY-noCells_0p0001
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p001.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p001
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p01.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p01
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p02.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p02
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p05.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p05
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p1.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p1
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p2.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p2
# mkdir ./output
# ./ovarian_follicle "./config/EGandGLY-withCells_0p5.xml" >> ./output/output_log
# mv ./output/ ./EGandGLY-noCells_0p5
# mkdir ./output
#
# ./ovarian_follicle "./config/GLY-noCells_0p0001.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p0001
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p001.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p001
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p01.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p01
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p02.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p02
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p05.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p05
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p1.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p1
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p2.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p2
# mkdir ./output
# ./ovarian_follicle "./config/GLY-noCells_0p5.xml" >> ./output/output_log 
# mv ./output/ ./GLY-withCells_0p5
# mkdir ./output


# ./ovarian_follicle "./config/HM05-withCells_0p0001.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p0001
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p001.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p001
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p01.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p01
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p02.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p02
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p05.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p05
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p1.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p1
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p2.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p2
# mkdir ./output
# ./ovarian_follicle "./config/HM05-withCells_0p5.xml" >> ./output/output_log 
# mv ./output/ ./HM05-noCells_0p5
# mkdir ./output
#
# ./ovarian_follicle "./config/HM1-withCells_0p0001.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p0001
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p001.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p001
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p01.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p01
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p02.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p02
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p05.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p05
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p1.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p1
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p2.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p2
# mkdir ./output
# ./ovarian_follicle "./config/HM1-withCells_0p5.xml" >> ./output/output_log 
# mv ./output/ ./HM1-noCells_0p5
# mkdir ./output
#
#
# ./ovarian_follicle "./config/HM2-withCells_0p0001.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p0001
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p001.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p001
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p01.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p01
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p02.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p02
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p05.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p05
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p1.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p1
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p2.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p2
# mkdir ./output
# ./ovarian_follicle "./config/HM2-withCells_0p5.xml" >> ./output/output_log 
# mv ./output/ ./HM2-noCells_0p5
# mkdir ./output
#
#
# ./ovarian_follicle "./config/HM5-withCells_0p0001.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p0001
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p001.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p001
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p01.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p01
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p02.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p02
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p05.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p05
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p1.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p1
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p2.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p2
# mkdir ./output
# ./ovarian_follicle "./config/HM5-withCells_0p5.xml" >> ./output/output_log 
# mv ./output/ ./HM5-noCells_0p5
# mkdir ./output
