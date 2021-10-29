#!/bin/bash
#
# Generates overview image of individual plants, segmented and original
#

IMAGE_DIR="processed"

for IMG_TYPE in seg image ; do
montage \
  -geometry 80x80+12+12 \
  -label '%t' \
  -pointsize 8 \
  $IMAGE_DIR/*/$IMG_TYPE/*.png \
  -frame 0 \
  -tile 12x \
  ${IMAGE_DIR}_${IMG_TYPE}_overview.png
done
