import cv2
import numpy as np
import os
import logging
import pdb
from skimage.util import img_as_float
from skimage import io
from skimage.measure import label, regionprops
from seg_slic import expand_mask
from PIL import Image

################################################################################
# Mask helper functions
################################################################################
def dilate_mask(ref_mask, dilate_size):
    # Prepare mask and image
    dilate_size = np.array(dilate_size, dtype=np.uint8)
    ref_mask = np.array(ref_mask, dtype=np.uint8)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (dilate_size, dilate_size))
    ref_mask_d = cv2.dilate(ref_mask, kernel)
    return ref_mask_d


def erode_mask(ref_mask, erode_size):
    erode_size = np.array(erode_size, dtype=np.uint8)
    ref_mask = np.array(ref_mask, dtype=np.uint8)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (erode_size, erode_size))
    ref_mask_d = cv2.erode(ref_mask, kernel)
    return ref_mask_d


def load_and_normalize_mask(mask_filename):
    """
    Load a mask file and normalize the indices
    Used by rapa_tray
    :param mask_filename: file name of the mask, should be an indexed png file
    :return: normalized mask, with the first index being 1, min and max indices
    """

    mask = np.array(Image.open(mask_filename))
    mask = np.array(mask, dtype='int16')
    current_min_idx = mask[np.where(mask > 0)].min()

    # Normalize mask
    mask = (mask + 1 - current_min_idx).copy()
    current_max_idx = mask.max()
    current_min_idx = mask[np.where(mask > 0)].min()

    # Check mask
    unique_values = np.unique(mask)
    position_indices = unique_values[np.where(unique_values > 0)]
    print(unique_values)
    is_mask_valid = np.array_equal(np.array([i for i in range(current_min_idx, current_max_idx + 1)]), position_indices)
    if not is_mask_valid:
        raise AssertionError('Mask error: position indices not valid. Please check continuity. Current set: {}'.format(position_indices))

    return mask, current_min_idx, current_max_idx


def do_grabcut(seg_image, img, gc_outer_iter, gc_mask_ext_size):
    ########################################################################
    # Refine segmentation
    # GrabCut postprocessing, from seg_slic.py
    if len(seg_image.shape) == 3:
        mask = (seg_image[:, : , 0] > 0 ) | (seg_image[:, : , 1] > 0 ) | (seg_image[:, : , 2] > 0 )
    else:
        mask = (seg_image > 0)
    image_for_grabcut = (img * 255).astype('uint8')
    if gc_outer_iter > 0:
        gc_mask = mask.astype('uint8')
        gc_mask[:,:] = cv2.GC_BGD
        gc_mask[np.where(mask == 0)] = cv2.GC_BGD
        gc_mask[np.where(mask == 1)] = cv2.GC_FGD
        gc_init = cv2.GC_INIT_WITH_MASK
        bgdmodel = np.zeros((1, 65), np.float64)
        fgdmodel = np.zeros((1, 65), np.float64)

        for iter in range(0, gc_outer_iter):
            print("GrabCut iteration " + str(iter))
            if iter == 0:
                gc_iter = 15
            else:
                gc_iter = 4
            gc_mask = expand_mask(gc_mask, gc_mask_ext_size)
            if sum(gc_mask.ravel()) > 100:
                # set_trace()
                cv2.grabCut(image_for_grabcut, gc_mask, None, bgdmodel, fgdmodel, gc_iter, gc_init)
            else:
                print("Mask size too small, no GrabCut!")
        mask2 = np.where((gc_mask == 2) | (gc_mask == 0), 0, 1).astype('uint8')
        seg_image = img * mask2[:, :, np.newaxis]
    return seg_image
    #
    ########################################################################




class Plate:
    """
    Contains steps for handling of the plate images
    """

    def __init__(self, plate_image=None, plant_mask=None, output_dir=None, plate_id=None,
                    area_of_first_mask=None):
        self.plate_image = None
        self.plant_mask = None
        self.output_dir = None if output_dir is None else output_dir
        self.plate_id = None if plate_id is None else plate_id
        self.area_of_mask = None if area_of_first_mask is None else area_of_first_mask
        self.plants = []  # Contains a dictionary with keys 'image', 'seg', 'green_px'
        self.logger = logging.getLogger('log')

        if plate_image is not None:
            self.load_image(plate_image)
        if plant_mask is not None:
            self.load_mask(plant_mask)

        if self.area_of_mask:
            self.pixels_to_area_factor = self.area_of_mask / \
                                        (self.plant_mask == 1).sum()
        else:
            self.pixels_to_area_factor = None


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

        self.plant_mask, _, _ = load_and_normalize_mask(mask_filename)


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

    def segment_image(self,input, output,
    L_range, a_range, b_range,
    kernel_size=8,
    remove_yellow_circles_at_border=False):
        img = io.imread(input, plugin='matplotlib')
        if output is None:
            output_filename = input_filename[:-4] + '_seg' + input_filename[-4:]
        else:
            output_filename = output
        img = img_as_float(img)
        img = img[:,:,0:3]

        from skimage.color import rgb2lab
        img_lab=rgb2lab(img)
        msk = (img_lab[:, :, 0] >= L_range[0]) * (img_lab[:, :, 0] < L_range[1]) *\
              (img_lab[:, :, 1] >= a_range[0]) * (img_lab[:, :, 1] < a_range[1]) *\
              (img_lab[:, :, 2] >= b_range[0]) * (img_lab[:, :, 2] < b_range[1])

        # clean up noise
        opened_msk = dilate_mask(erode_mask(msk, kernel_size), kernel_size)

        seg_image = img * opened_msk[:, :, np.newaxis]

        ########################################################################
        # Remove yellow circles at borders specific to Talias pots
        # Parameters:
        from skimage.morphology import disk
        from skimage.morphology import binary_dilation
        from skimage.measure import label,regionprops
        fill_holes_kernel_size = 20 # For filling holes inside the regions
        corner_dist = 100 # removal distance threshold to one of the corners
        overlap_fraction = 0.75
        if remove_yellow_circles_at_border:
            # Define color ranges, determined by threshold analyzer
            yellow_circle_b_range = (60, 127)

            # Create mask
            cmask = (img_lab[:, :, 2] >= yellow_circle_b_range[0]) *\
                    (img_lab[:, :, 2] < yellow_circle_b_range[1])

            # Fill up holes
            # Kernel size of 20 is roughly the width of the pen strocke on the
            # circular label
            dilated_cmask = dilate_mask(cmask, fill_holes_kernel_size)

            # Removal decision: if the region is within corner_dist px distance
            #to one of the corners, remove. corner_dist is set above.
            removal_zones_mask = dilated_cmask * 0
            removal_zones_mask[0:corner_dist, 0:corner_dist] = 1
            removal_zones_mask[0:corner_dist, -corner_dist:] = 1
            removal_zones_mask[-corner_dist:, 0:corner_dist] = 1
            removal_zones_mask[-corner_dist:, -corner_dist:] = 1

            # Process regions
            labeled_regions = label(dilated_cmask)
            for region in regionprops(labeled_regions):
                # Is the region completely inside the removal zones?
                current_region_msk = (labeled_regions == region.label)
                if np.sum(current_region_msk) >= np.sum(current_region_msk * removal_zones_mask) * overlap_fraction:
                   seg_image = (current_region_msk == 0)[:,:,np.newaxis] * \
                   seg_image[:, :, :]
        #
        ########################################################################

        ########################################################################
        # Refine segmentation
        gc_outer_iter = 50
        gc_mask_ext_size = 15
        # GrabCut postprocessing, from seg_slic.py
        seg_image = do_grabcut(seg_image, img, gc_outer_iter, gc_mask_ext_size)
        #
        ########################################################################

        ########################################################################
        # Add purple leaves
        # postprocess_purple_radius = int(img.shape[0]/4)
        # purple_leaves_b_range = (-127, -10)
        # purple_leaves_mask = (img_lab[:, :, 2] >= purple_leaves_b_range[0]) *\
        #                     (img_lab[:, :, 2] < purple_leaves_b_range[1])
        # circular_mask = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (postprocess_purple_radius * 2,
        #                                                               postprocess_purple_radius * 2))
        # pad_widths = [int((img.shape[i] - postprocess_purple_radius * 2)/2) + 1 for i in [0, 1]]
        # circular_mask = np.pad(circular_mask, pad_width=((pad_widths[0], pad_widths[0]), (pad_widths[1], pad_widths[1])),
        #        mode='constant', constant_values=0)
        # circular_mask = circular_mask[0:img.shape[0], 0:img.shape[1]]
        # labeled_regions = label(purple_leaves_mask)
        # for region in regionprops(labeled_regions):
        #     current_region_msk = (labeled_regions == region.label)
        #     if np.sum(current_region_msk) != np.sum(current_region_msk * circular_mask):
        #         purple_leaves_mask = (1 - current_region_msk) * purple_leaves_mask
        # purple_leaves_seg_image = do_grabcut(purple_leaves_mask, img, gc_outer_iter, gc_mask_ext_size)
        # purple_leaves_mask_post_gc = (purple_leaves_seg_image[:, : , 0] > 0 ) | (purple_leaves_seg_image[:, : , 1] > 0 ) | (purple_leaves_seg_image[:, : , 2] > 0)
        # mask_to_three_channels = lambda x: np.stack((x,)*3, axis=-1)
        # purple_indices = np.where(mask_to_three_channels(purple_leaves_mask_post_gc) > 0)
        # seg_image[purple_indices] = img[purple_indices]
        ########################################################################
        # Clean up GrabCut residue, remove unsaturated pixels
        unsaturated_mask = (img_lab[:, :, 0] >= 20) * (img_lab[:, :, 0] <= 100) *\
                           (img_lab[:, :, 1] >= -10) * (img_lab[:, :, 1] < 30) *\
                           (img_lab[:, :, 2] >= -10) * (img_lab[:, :, 2] < 30)
        seg_image = (1 - unsaturated_mask)[:,:,np.newaxis] * seg_image
        #
        ########################################################################

        io.imsave(output, seg_image)
        return seg_image




    def segment_plants(self):
        for i,p in enumerate(self.plants):
            target_dir = self.get_plant_filename(i, 'seg')[1]
            if not os.path.isdir(target_dir):
                os.makedirs(target_dir)
            target_file = self.get_plant_filename(i, 'seg')[0]
            image_file = self.get_plant_filename(i, 'image')[0]
            segmented_image = self.segment_image(input=image_file, output=target_file,
            L_range=(-127,127),
            a_range=(-127,-10),
            b_range=(-127,127),
            kernel_size=8,
            remove_yellow_circles_at_border=True)
            p['seg'] = cv2.imread(target_file)


    def calc_area(self):
        for i,p in enumerate(self.plants):
            p["area"] = ((p["seg"][:,:,0]>0) * (p["seg"][:,:,1] > 0) * (p["seg"][:,:,2] > 0)).sum() 

            
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
                    if self.pixels_to_area_factor is None:
                        of.write('{0},{1:03},{2}\n'.format(self.plate_id, i, p['area']))
                    else:
                        of.write('{0},{1:03},{2},{3}\n'.format(self.plate_id, i, p['area'], self.pixels_to_area_factor * p['area']))

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
