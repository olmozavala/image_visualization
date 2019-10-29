import SimpleITK as sitk
from os.path import join
from img_viz.constants import SliceMode, PlaneTypes

from img_viz.medical import MedicalImageVisualizer
from img_viz.maskrcnn import MaskRCNNVisualizer


def main():
    input_folder = join('..', 'test_data', 'input_medical', 'Case-0001')
    output_folder = join('..', 'test_data', 'output_medical')
    itk_img = sitk.ReadImage(join(input_folder,'img_tra.nrrd'))
    ctr1 = sitk.ReadImage(join(input_folder,'ctr_Skin.nrrd'))
    ctr2 = sitk.ReadImage(join(input_folder,'ctr_KidneyBilat.nrrd'))

    viz_obj = MedicalImageVisualizer(output_folder=output_folder, disp_images=False)
    # ================= Single Slice (middle) for all axes from ITK images ==============
    print("Making plots for a single slice...")
    title='Single Slice'
    viz_obj.plot_img_and_ctrs_itk(itk_img=itk_img, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
                                  title=title, plane=PlaneTypes.ALL, draw_only_ctrs=True,
                                  file_name_prefix='AllPlanes', labels=['Prostate', 'PZ'])
    viz_obj.plot_img_and_ctrs_itk(itk_img=itk_img, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
                                  title=title, plane=PlaneTypes.AXIAL, draw_only_ctrs=True,
                                  file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])
    viz_obj.plot_img_and_ctrs_itk(itk_img=itk_img, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
                                  title=title, plane=PlaneTypes.CORONAL, draw_only_ctrs=True,
                                  file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])
    viz_obj.plot_img_and_ctrs_itk(itk_img=itk_img, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
                                  title=title, plane=PlaneTypes.SAGITTAL, draw_only_ctrs=True,
                                  file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])

    # ================= Multiple Slices  for single axe from ITK images ==============
    print("Making plots for a multiple slices...")
    slices = range(0,100, 10)
    title = 'Multiple Slices only contours'
    viz_obj.plot_img_and_ctrs_itk(itk_img=itk_img, itk_ctrs=[ctr1, ctr2], slices=slices,
                                  title=title, plane=PlaneTypes.CORONAL, draw_only_ctrs=True,
                                  file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])
    title = 'Multiple Slices all'
    viz_obj.plot_img_and_ctrs_itk(itk_img=itk_img, itk_ctrs=[ctr1, ctr2], slices=slices,
                                  title=title, plane=PlaneTypes.SAGITTAL, draw_only_ctrs=False,
                                  file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])

    # ================= Plot locations ===============
    print("Making plots for a 3D locations...")
    loc = [[40,100,200], [20,130,208]]
    viz_obj.plot_locations_np(c_img=sitk.GetArrayFromImage(itk_img), locations=loc,
                              file_name_prefix='SinglePlane')


    # ================= Plot Statistics ===============
    print("Making plots for a statistics...")
    data = [{'c1': 10, 'c2': 14, 'c3': 10},
            {'c1': 12, 'c2': 8, 'c3': 13},
            {'c1': 11, 'c2': 9, 'c3': 10}]
    legends = ['First', 'Second', 'Third']
    viz_obj.plot_multiple_bar_plots(data, file_name='Bar_plots.png', title='Cool plot', legends=legends)

    viz_obj.plot_histograms_from_imgs([itk_img, itk_img],file_name='img_histo.jpg',titles=['Normal', 'Normalized'])


if __name__ == '__main__':
    main()
