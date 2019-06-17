#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import logging, coloredlogs
#import plate
sys.path.append("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/code")
import plate

if __name__ == "__main__":

    ############################ Parameter ############################
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input_image', required=True,
                    default=None, help='Input image file')
    ap.add_argument('-pid', '--plate_id', required=True,
                    default=None, help='Plate ID, used as prefix for file numbering')
    ap.add_argument('-o', '--output_dir', required=False,
                    default='.', help='Output directory for the plate subfolder, default is .')
    ap.add_argument('-m', '--mask', required=True,
                    default=None, help='Mask for plants, should be an indexed png file with ascending indices.')
    args = vars(ap.parse_args())
    logger = logging.getLogger('log')
    if not os.path.isdir(args['output_dir']):
        try:
            os.makedirs(args['output_dir'])
        except OSError:
            logger.critical('Could not create directory {}!'.format(args['block_output_dir']))
            sys.exit(1)
    if not os.path.isfile(args['input_image']):
        logger.critical('Input image {} does not exist!'.format(args['input_image']))
        sys.exit(1)


    ############################ Processing #############################
    #pl = plate.Plate(plate_image='./input_images/plate10.jpg', plant_mask='./input_images/mask.png',output_dir='./processed_images', plate_id='plate10')

    pl = plate.Plate(plate_image=args['input_image'], plant_mask=args['mask'],
                     output_dir=args['output_dir'], plate_id=args['plate_id'])
    pl.extract_plants_from_images()
    pl.save_everything()
    pl.segment_plants()
    pl.calc_area()
    pl.save_everything()





