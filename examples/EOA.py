import os
import xarray as xr
os.environ['PROJ_LIB'] = '/home/olmozavala/anaconda3/envs/eoas/share/proj'

from os.path import join
from img_viz.constants import SliceMode, PlaneTypes

from img_viz.eoa_viz import EOAImageVisualizer

from netCDF4 import Dataset


def main():
    input_folder = join('..', 'test_data', 'input_eoas')
    output_folder = join('..', 'test_data', 'output_eoas')

    nc_ds_u = Dataset(join(input_folder,'ex1_u.nc'), "r+", format="NETCDF4")
    nc_ds_v = Dataset(join(input_folder,'ex1_v.nc'), "r+", format="NETCDF4")

    xa_ds_u = xr.open_dataset(join(input_folder,'ex1_u.nc'))

    viz_obj = EOAImageVisualizer(output_folder=output_folder, disp_images=True)

    # ================= SUMMARY of file ====================
    viz_obj.nc_summary(nc_ds_u)
    viz_obj.nc_summary(nc_ds_v)

    # ================= 4D data ============
    # --------- NetCDF Maps  -----
    var_names = ['u']
    timesteps = [0]
    z_axis_levels = [0,1]
    var_names = ['v']
    # viz_obj.plot_4d_data_ncds_map(nc_ds_v, var_names=var_names, z_levels=z_axis_levels, timesteps=timesteps,
    #                           title='Title V', file_name_prefix='Test_V')

    # --------- Xarray  Maps ---
    var_names = ['u']
    viz_obj.plot_4d_data_xarray_map(xa_ds_u, var_names=var_names, z_levels=z_axis_levels, timesteps=timesteps,
                                    title='Title U', file_name_prefix='Xarray_map_U')

    # --------- NetCDF (Normally dont use this one) -----
    # viz_obj.plot_4d_data_ncds(nc_ds_u, var_names=var_names, z_levels=z_axis_levels, timesteps=timesteps,
    #                           title='Title U', file_name_prefix='Test_U')
    # var_names = ['v']
    # viz_obj.plot_4d_data_ncds(nc_ds_v, var_names=var_names, z_levels=z_axis_levels, timesteps=timesteps,
    #                           title='Title V', file_name_prefix='Test_V')

    # ================= 3D data ============
    # cur_var = nc_ds_u.variables['u'][0]  # Getting only the first time
    # viz_obj.plot_3d_data_singlevar_np(cur_var, z_axis_levels, title='Numpy3D',
    #                                   file_name_prefix='numpy_3d')

if __name__ == '__main__':
    main()

##

