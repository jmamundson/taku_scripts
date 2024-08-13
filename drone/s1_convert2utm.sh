#/bin/bash
#
# convert lat/lon to UTM Zone 8N; use -tr to specify pixel size

gdalwarp -t_srs epsg:32608 -tr 0.2 0.2 $1 $2
