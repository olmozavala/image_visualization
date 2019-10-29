import SimpleITK as sitk
from os.path import join
from img_viz.constants import SliceMode, PlaneTypes

from img_viz.medical import MedicalImageVisualizer
from img_viz.maskrcnn import MaskRCNNVisualizer


def main():
    # ================= Generate Images for MaskRCNN ===============
    input_folder = join('..','test_data', 'input_medical')
    output_folder = join('..', 'test_data', 'output_medical')
    viz_maskrcnn_obj = MaskRCNNVisualizer(output_folder=output_folder,
                                          input_folder=input_folder,
                                          disp_images=False)

    viz_maskrcnn_obj.plot_for_mask_r_cnn_input(img_name='img_tra.nrrd',
                                               ctr_names=['ctr_Pancreas.nrrd', 'ctr_Skin.nrrd', 'ctr_Stomach.nrrd'],
                                               slices=SliceMode.MIDDLE, plane=PlaneTypes.AXIAL, plot_empty_ctr=True)
    viz_maskrcnn_obj.plot_for_mask_r_cnn_input(img_name='img_tra.nrrd',
                                               ctr_names=['ctr_Pancreas.nrrd', 'ctr_Skin.nrrd', 'ctr_Stomach.nrrd'],
                                               slices=SliceMode.MIDDLE, plane=PlaneTypes.SAGITTAL, plot_empty_ctr=True)

if __name__ == '__main__':
    main()
