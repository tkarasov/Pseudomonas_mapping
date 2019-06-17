#!/bin/bash


#ssh -XC tkarasov@burrito.eb.local
bash
source activate cv
source /ebio/abt6_rapa/software/snapshots/rapa_flir/src/pipeline/rapa_env.bash
cd /ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/
declare -A IDS
for i in {1..32};
do
IDS[plate$i.jpg]=plate$i
done
#IDS["IMG_8985.JPG"]=23
#IDS["IMG_8987.JPG"]=24
#IDS["IMG_8994.JPG"]=27
MASK=mask_day7.png
OUTPUT_DIR=processed_images

for IMG in ${!IDS[@]}; do
    ./code/process_plate.py --input_image ./input_images/$IMG --plate_id ${IDS[${IMG}]} -m ./input_images/$MASK -o processed_images
done


#now map the pixels to the actual strain
python /ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/code/collate_plates.py

#this output a file 'cell_pixel_strain.txt'
