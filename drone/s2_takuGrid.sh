#!/bin/bash
#
# create same grid for all products

base=/hdd2/taku/uav/surveys

#for file in $(ls $base'/01_utm/*_2025*.tif'); do

output=$base'/02_takuGrid/'$(basename $1)

xmin=554000 
xmax=559000
ymin=6473000
ymax=6477000

gdalwarp -te $xmin $ymin $xmax $ymax -tr 1  1 $1 $output;

#done
