#!/bin/bash
# Version: 20190725_aduque

##############################################################################
# User parameters, adjust as needed
##############################################################################
IMAGE_DIR=$(pwd) # Adjust if needed
MASK="$IMAGE_DIR/mask.png" # Mask, adjust if necessary
OUTPUT_DIR="$(pwd)/processed"
NUMBER_OF_JOBS=6  # Don't be too greedy

##############################################################################
# Internal variables
##############################################################################
CURRENT_DIR=$(pwd)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROCESS_PLATE=${SCRIPT_DIR}/process_plate


##############################################################################
# Processing
##############################################################################
cd $IMAGE_DIR
echo "Running segmentations"

# Use GNU parallel for parallel processing of the images
# {} is the filename
# plate id is set to {.}, which means the basename without the extension:
# input/test_plate.JPG becomes test_plate
parallel --nn --ungroup -j ${NUMBER_OF_JOBS} \
   "$PROCESS_PLATE \
     --input_image {} \
      --plate_id {.} \
      --output_dir $OUTPUT_DIR \
      --mask $MASK" ::: *.JPG
 
cd $CURRENT_DIR

echo "Creating csv and image overview"
${SCRIPT_DIR}/create_area_csv_file.bash
${SCRIPT_DIR}/create_overview_image.bash
