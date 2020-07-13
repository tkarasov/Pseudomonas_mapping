import rapa_utils as ru
import numpy as np
import os
import pdb
import sys
import cv2
import logging, coloredlogs
class Plate:
    """
    Contains steps for handling of the plate images
    """

    def __init__(self, plate_image=None, plant_mask=None, output_dir=None, plate_id=None):
        self.plate_image = None
        self.plant_mask = None
        self.output_dir = None if output_dir is None else output_dir
        self.plate_id = None if plate_id is None else plate_id
        self.plants = []  # Contains a dictionary with keys 'image', 'seg', 'green_px'
        self.logger = logging.getLogger('log')

        if plate_image is not None:
            self.load_image(plate_image)
        if plant_mask is not None:
            self.load_mask(plant_mask)


    def load_image(self, image_filename):
        """
        Loads a plate image
        :param image_filename: image file name
        :return:
        """
        self.plate_image = cv2.imread(image_filename)


    def load_mask(self, mask_filename):
        """
        Load a mask image, which is an indexed png file where the masks for each
        plant are idenfied by ascending indices.
        :param mask_filename: mask file name
        :return:
        """

        self.plant_mask, _, _ = ru.load_and_normalize_mask(mask_filename)


    def extract_plants_from_images(self):
        """
        Extracts the plants from the loaded image
        :return:
        """

        assert self.plate_image is not None
        assert self.plant_mask is not None
        for p in range(1, np.max(self.plant_mask[:])+1):
            curr_plant = {}
            curr_plant_msk = self.plant_mask == p
            p_img = self.plate_image * np.uint8(curr_plant_msk[:, :, np.newaxis])
            # get bounding box
            m_ind = np.where(curr_plant_msk == True)
            curr_plant['image'] = (self.plate_image[m_ind[0].min():m_ind[0].max(),
                                 m_ind[1].min():m_ind[1].max(),:])
            self.plants.append(curr_plant)


    def segment_plants(self):
        import seg_slic
        for i,p in enumerate(self.plants):
            target_dir = self.get_plant_filename(i, 'seg')[1]
            if not os.path.isdir(target_dir):
                os.makedirs(target_dir)
            target_file = self.get_plant_filename(i, 'seg')[0]
            image_file = self.get_plant_filename(i, 'image')[0]
            seg_slic.segment_image(image_file, 0, 0, 10, 'rgb', target_file,
                                   None, None, None, None, 'default')
            p['seg'] = cv2.imread(target_file)

    def calc_area(self):
        for i,p in enumerate(self.plants):
            p['area'] = (p['seg']>0).sum()

    def save_everything(self):
        """
        Saves all currently present data
        :param output_dir: Target directory
        :return:
        """

        assert os.path.isdir(self.output_dir)
        assert self.plate_id is not None

        plate_dir = os.path.join(self.output_dir, self.plate_id)
        if not os.path.isdir(plate_dir):
            os.makedirs(plate_dir)

        for i,p in enumerate(self.plants):
            for k in ['image', 'seg']:
                if k in p:
                    imdir = os.path.join(plate_dir, k)
                    if not os.path.isdir(imdir):
                        os.makedirs(imdir)
                    filename, _ = self.get_plant_filename(i, k)
                    cv2.imwrite(filename, p[k])
        if 'area' in self.plants[0]:
            with open(os.path.join(imdir,'..', 'area.txt'), 'w') as of:
                for i,p in enumerate(self.plants):
                    of.write('{0},{1:03},{2}\n'.format(self.plate_id, i, p['area']))


    def get_plant_filename(self, plant_number, key):
        """
        return the plant filename
        :param plant_number:
        :param key:
        :return:
        """
        plate_dir = os.path.join(self.output_dir, self.plate_id)
        imdir = os.path.join(plate_dir, key)
        filename = '{0}_{1:03d}_{2}.png'.format(self.plate_id, plant_number, key)
        return os.path.join(imdir,filename),imdir


    def write_overview_image(self):
        self.logger.critical('Writes an overview image with area underneath. not impemented!')
        raise NotImplementedError




