#!/bin/bash
#
# convert DEMs into hillshades for feature tracking
base=/hdd/taku/uav/surveys

for file in $(ls $base'/02_takuGrid/taku_DEM'*.tif); do

output=$base'/03_hillshades/'$(basename $file)    

gdaldem hillshade $file $output;

done
