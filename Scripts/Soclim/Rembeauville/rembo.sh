#!/bin/sh

floatID=$1
floatWMO=$2
dac=$3

echo $floatID
echo $floatWMO
echo $dac

R $floatID $floatWMO $dac --vanilla  < WMOtoRAMBOVILLE_format.R 
python prediction_AP_VER4.py $floatID
R $floatID $floatWMO $dac --vanilla < RAMBOVILLE_out_plot.R



#./rembo.sh 049b 6901585 coriolis
#./rembo.sh 036b 6901583 coriolis
#./rembo.sh 037c 6901004 coriolis
#./rembo.sh 107c 6902739 coriolis
#./rembo.sh 104c 6902738 coriolis
