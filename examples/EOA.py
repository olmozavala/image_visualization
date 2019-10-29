from os.path import join
from img_viz.constants import SliceMode, PlaneTypes

from img_viz.eoa_viz import EOAImageVisualizer

from netCDF4 import Dataset

def main():
    input_folder = join('..', 'test_data', 'input_eoas')
    output_folder = join('..', 'test_data', 'output_eoas')

    nc_ds_u = Dataset(join(input_folder,'ex1_u.nc'), "r+", format="NETCDF4")
    nc_ds_v = Dataset(join(input_folder,'ex1_v.nc'), "r+", format="NETCDF4")

    viz_obj = EOAImageVisualizer(output_folder=output_folder, disp_images=True)

    # ================= Summary of the file ====================
    viz_obj.nc_summary(nc_ds_u)

    # ================= Plot single frame ============
    var_names = ['u']
    timesteps = [0]
    z_axis_levels= [0,1,2]
    viz_obj.plot_4d_data_np(nc_ds_u, var_names, timesteps, z_axis_levels, 'Example Title', 'Test_U')
    # ================= Single Slice (middle) for all axes from ITK images ==============
    # print("Making plots for a single slice...")
    # title='Single Slice'
    # viz_obj.plot_img_and_ctrs_itk(nc_ds=nc_ds, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
    #                               title=title, plane=PlaneTypes.ALL, draw_only_ctrs=True,
    #                               file_name_prefix='AllPlanes', labels=['Prostate', 'PZ'])
    # viz_obj.plot_img_and_ctrs_itk(nc_ds=nc_ds, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
    #                               title=title, plane=PlaneTypes.AXIAL, draw_only_ctrs=True,
    #                               file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])
    # viz_obj.plot_img_and_ctrs_itk(nc_ds=nc_ds, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
    #                               title=title, plane=PlaneTypes.CORONAL, draw_only_ctrs=True,
    #                               file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])
    # viz_obj.plot_img_and_ctrs_itk(nc_ds=nc_ds, itk_ctrs=[ctr1, ctr2], slices=SliceMode.MIDDLE,
    #                               title=title, plane=PlaneTypes.SAGITTAL, draw_only_ctrs=True,
    #                               file_name_prefix='SinglePlane', labels=['Prostate', 'PZ'])

if __name__ == '__main__':
    main()
