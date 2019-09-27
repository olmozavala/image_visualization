import numpy as np
import datetime
import os
from os.path import join
import SimpleITK as sitk
from img_viz.medical import MedicalImageVisualizer
from img_viz.constants import SliceMode, PlaneTypes


# This code can be used as an example when we want to generate images from a preproc folder.
def main():

    input_folder = '/media/osz1/DATA/DATA/GE/Preproc'
    output_folder = '/media/osz1/DATA/DATA/IMAGES/TEST'

    ctr_names = ['roi_ctr_pro.nrrd', 'roi_ctr_pz.nrrd']
    img_names = ['roi_tra.nrrd', 'roi_sag.nrrd']
    labels = ['Prostate', 'PZ']

    plot_multiple_cases(input_folder, output_folder, img_names, ctr_names, labels)


def plot_multiple_cases(input_folder, output_folder, img_names, ctr_names, labels):
    """
    It simply calls the plot function for every image and case found in the input folder
    :param input_folder:
    :param output_folder:
    :param img_names:
    :param ctr_names:
    :return:
    """

    # Initialize Object
    viz_obj = MedicalImageVisualizer(output_folder=output_folder, disp_images=False)

    all_cases = os.listdir(input_folder)
    # Iterate over each image
    # Iterate over each case
    for c_case in all_cases:
        try:
            print(F'***** {c_case} ********')
            # Create image name
            out_img_file_name = F'{c_case}'
            ctrs_itk = []
            imgs_itk = []
            for c_img_name in img_names:
                img_itk = sitk.ReadImage(join(input_folder, c_case, c_img_name))
                imgs_itk.append(img_itk)

            # Read all the desired contours
            for c_ctr_name in ctr_names:
                ctr_itk = sitk.ReadImage(join(input_folder, c_case, c_ctr_name))
                ctrs_itk.append(ctr_itk)

            # Plot
            viz_obj.plot_imgs_and_ctrs_itk(imgs_itk, ctrs_itk, slices=SliceMode.MIDDLE,
                                           draw_only_ctrs=False, file_name_prefix=out_img_file_name,
                                           labels=labels)

        except Exception as e:
            print(F' Failed for case {c_case}, error: {e}')


if __name__ == '__main__':
    main()
