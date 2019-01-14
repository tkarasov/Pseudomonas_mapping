#!/usr/bin/env python3

#this script takes the output from the process_plate.py pipeline, and assigns strain identities to each image.
import os
import sys
import glob



if __name__ == "__main__":

    ############################ Parameter ############################
    '''
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--image_directory', required=True,
                    default=None, help='Input image file')
    ap.add_argument('-o', '--output_dir', required=False,
                    default='.', help='Output directory for concatenated data, default is .')
    ap.add_argument('-m', '--metadata', required=False,
                    default='.', help='Meta data for the plate layout, default is .')

    args = vars(ap.parse_args())
    logger = logging.getLogger('log')
    if not os.path.isdir(args['output_dir']):
        try:
            os.makedirs(args['output_dir'])
        except OSError:
            logger.critical('Could not create directory {}!'.format(args['block_output_dir']))
            sys.exit(1)
    if not os.path.isdir(args['output_dir']):
        logger.critical('Input image {} does not exist!'.format(args['input_image']))
        sys.exit(1)
'''
    ############################ Processing #############################
    meta = [line.strip().split() for line in open("/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/image_analysis/plate_position_orig_strain_map_exact.txt").readlines()]
    image_directory = '/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/processed_images'
    
    #build dictionary with every plate location and its strain
    plate_dict = {}
    for line in meta:
    	plate_dict[(line[0], line[1])]=[line[3]]

    #build dictionary of plate location mapping
    plate_num = {}
    plate_num['000'] = "A1"
    plate_num['001'] = "A2"
    plate_num['002'] = "A3"
    plate_num['003'] = "A4"
    plate_num['004'] = "A5"
    plate_num['005'] = "A6"
    plate_num['006'] = "B1"
    plate_num['007'] = "B2"
    plate_num['008'] = "B3"
    plate_num['009'] = "B4"
    plate_num['010'] = "B5"
    plate_num['011'] = "B6"
    plate_num['012'] = "C1"
    plate_num['013'] = "C2"
    plate_num['014'] = "C3"
    plate_num['015'] = "C4"
    plate_num['016'] = "C5"
    plate_num['017'] = "C6"
    plate_num['018'] = "D1"
    plate_num['019'] = "D2"
    plate_num['020'] = "D3"
    plate_num['021'] = "D4"
    plate_num['022'] = "D5"
    plate_num['023'] = "D6"

    #read in every read file from image_directory, put number of pixels into plate_dict
    for dir in glob.iglob(image_directory+"/plate*"):
    	area = [line.strip().split(',') for line in open(dir+"/area.txt").readlines()]
    	for cell in area:
    		position = plate_num[str(cell[1])]
    		plate = cell[0]
    		area_count = cell[2]
    		plate_dict[(plate, position)].append(cell[2])


    with open('/ebio/abt6_projects9/Pseudomonas_diversity/Pseudomonas_mapping/data/infection_experiments/june_2018_day7/image_analysis/cell_pixel_strain.txt', 'w') as f:
    	for key, value in plate_dict.items():
    		f.write('{} {}\n'.format(' '.join(key), ' '.join(value)))

