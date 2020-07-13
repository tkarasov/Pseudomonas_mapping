#!/bin/bash
# 
# Generates csv file for pixel counts of all plates
# 
IMAGE_DIR="processed"

CSV_FILE=${IMAGE_DIR}_area_in_pixels.csv

echo "Plate,Pot,Green_pixels" > $CSV_FILE
cat $IMAGE_DIR/*/area.txt  >> $CSV_FILE

