import SimpleITK as sitk
from os.path import join
from scipy import ndimage

from img_viz.viz3d import ImageVisualizer3D


def main():
    input_folder = join('..', 'test_data', 'input', 'Case-0001')
    output_folder = join('..', 'test_data', 'output')
    ctr_itk = sitk.ReadImage(join(input_folder, 'ctr_KidneyBilat.nrrd'))

    viz_obj = ImageVisualizer3D(output_folder=output_folder, open_browser=False)
    print("Making surface....")
    viz_obj.plot_surface_itk(ctr_itk, title='Surface Title', file_name='Surface_3D_example')
    print("Done!!!")

    print("Making scatter 3D....")
    # In this case we want only the 'contour' or it will take a looooooooooooooong time
    ctr_np = sitk.GetArrayFromImage(ctr_itk)
    ctr_eroded_np = ndimage.morphology.binary_erosion(ctr_np)
    ctr_boundary_np = ctr_np - ctr_eroded_np
    viz_obj.plot_scatter_np(ctr_boundary_np, title='Scatter Title', file_name='Scatter_3D_example')
    print("Done!!!")


if __name__ == '__main__':
    main()
