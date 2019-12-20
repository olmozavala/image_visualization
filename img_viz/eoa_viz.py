from os import listdir
from os.path import join

import matplotlib.pyplot as plt
from img_viz.common import *
from img_viz.constants import SliceMode, PlaneTypes
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

class EOAImageVisualizer:
    """This class makes plenty of plots from netcdf datasets or numpy arrays """
    _disp_images = True
    _output_folder = 'output_medical'
    _COLORS = ['y', 'r', 'c', 'b', 'g', 'w', 'k', 'y', 'r', 'c', 'b', 'g', 'w', 'k']

    def __init__(self, **kwargs):
        # All the arguments that are passed to the constructor of the class MUST have its name on it.
        for arg_name, arg_value in kwargs.items():
            self.__dict__["_" + arg_name] = arg_value

    def __getattr__(self, attr):
        '''Generic getter for all the properties of the class'''
        return self.__dict__["_" + attr]

    def __setattr__(self, attr, value):
        '''Generic setter for all the properties of the class'''
        self.__dict__["_" + attr] = value

    def _close_figure(self):
        """Depending on what is disp_images, the figures are displayed or just closed"""
        if self.disp_images:
            plt.show()
        else:
            plt.close()

    def getExtent(self, lats, lons):
        minLat = np.amin(lats)
        maxLat = np.amax(lats)
        minLon = np.amin(lons)
        maxLon = np.amax(lons)
        bbox = (minLon, maxLon, minLat, maxLat)
        return bbox

    def xr_summary(self, ds):
        """ Prints a summary of the netcdf (global attributes, variables, etc)
        :param ds:
        :return:
        """
        print("\n========== Global attributes =========")
        for name in ds.attrs:
            print(F"{name} = {getattr(ds, name)}")

        print("\n========== Dimensions =========")
        for name in ds.dims:
            print(F"{name}: {ds[name].shape}")

        print("\n========== Coordinates =========")
        for name in ds.coords:
            print(F"{name}: {ds[name].shape}")

        print("\n========== Variables =========")
        for cur_variable_name in ds.variables:
            cur_var = ds[cur_variable_name]
            print(F"{cur_variable_name}: {cur_var.dims} {cur_var.shape}")

    def nc_summary(self, ds):
        """ Prints a summary of the netcdf (global attributes, variables, etc)
        :param ds: 
        :return: 
        """
        print("\n========== Global attributes =========")
        for name in ds.ncattrs():
            print(F"{name} = {getattr(ds, name)}")

        print("\n========== Variables =========")
        netCDFvars = ds.variables
        for cur_variable_name in netCDFvars.keys():
            cur_var = ds.variables[cur_variable_name]
            print(F"Dimensions for {cur_variable_name}: {cur_var.dimensions} {cur_var.shape}")

    def add_roads(self, ax):
        # Names come from: https://www.naturalearthdata.com/features/
        # -- Add states
        roads = cfeature.NaturalEarthFeature(
            category='cultural',
            name='roads',
            scale='10m',
            facecolor='none')

        ax.add_feature(roads, edgecolor='black')
        return ax

    def add_states(self, ax):
        # Names come from: https://www.naturalearthdata.com/features/
        # -- Add states
        states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

        ax.add_feature(states_provinces, edgecolor='gray')
        return ax
    # ====================================== 4D Data =====================================

    def plot_4d_data_xarray_map(self, ds, var_names: list, z_levels: list, timesteps: list, title='',
                          file_name_prefix='', cmap='viridis', proj=ccrs.PlateCarree(), lonvar='Longitude', latvar='Latitude',
                                z_levelvar='Depth'):
        """
        Plots multiple z_levels from a NetCDF file (4D data). It plots the results in a map.
        http://xarray.pydata.org/en/stable/plotting.html#maps
        """
        projection = proj
        lon = ds[lonvar]
        lat = ds[latvar]

        checkFolder(self._output_folder)
        for c_slice in z_levels:
            for c_time in timesteps:
                plt.subplots(1, len(var_names), squeeze=True, figsize=(8*len(var_names), 8))
                for idx_var, c_var_name in enumerate(var_names):
                    # Depth selection, not sure if the name is always the same
                    c_depth = ds[z_levelvar].values[c_slice]
                    # Interpolates data to the nearest requested depth
                    cur_var = ds[c_var_name].sel(**{z_levelvar: c_depth}, method='nearest')
                    ax = plt.subplot(1, len(var_names), idx_var+1, projection=projection)
                    # ------------------ MAP STUFF -------------------------
                    # https://rabernat.github.io/research_computing_2018/maps-with-cartopy.html
                    cur_var.plot(ax=ax, transform=projection, cbar_kwargs={'shrink': 0.4})
                    # ax.stock_img()  # Draws a basic topography of the world
                    ax.coastlines(resolution='50m')  # Draws the coastline
                    # ax = self.add_states(ax)
                    ax = self.add_roads(ax)
                    ax.add_feature(cartopy.feature.OCEAN)
                    ax.add_feature(cartopy.feature.LAND, edgecolor='black')
                    ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
                    ax.add_feature(cartopy.feature.RIVERS)
                    ax.gridlines()
                    # ------------------ MAP STUFF -------------------------
                    c_title = F'{c_var_name} {title} Z-level:{c_slice} Time:{c_time} '
                    plt.title(c_title, fontsize=20)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()

    def plot_4d_data_ncds_map(self, ds, var_names: list, z_levels: list, timesteps: list, title='',
                          file_name_prefix='', cmap='viridis', proj=ccrs.PlateCarree(), lonvar='Longitude', latvar='Latitude'):
        """
        Plots multiple z_levels from a NetCDF file (4D data). It plots the results in a map.
        """
        lon = ds.variables[lonvar][:]
        lat = ds.variables[latvar][:]
        projection = proj

        checkFolder(self._output_folder)
        for c_slice in z_levels:
            for c_time in timesteps:
                plt.subplots(1, len(var_names), squeeze=True, figsize=(8*len(var_names), 8))
                for idx_var, c_var_name in enumerate(var_names):
                    cur_var = ds.variables.get(c_var_name)
                    ax = plt.subplot(1, len(var_names), idx_var+1, projection=projection)
                    # ------------------ MAP STUFF -------------------------
                    # https://rabernat.github.io/research_computing_2018/maps-with-cartopy.html
                    ax.stock_img()  # Draws a basic topography
                    ax.coastlines(resolution='50m')  # Draws the coastline
                    ax.add_feature(cartopy.feature.OCEAN)
                    ax.add_feature(cartopy.feature.LAND, edgecolor='black')
                    ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
                    ax.add_feature(cartopy.feature.RIVERS)
                    ax.gridlines()
                    # ------------------ MAP STUFF -------------------------
                    # plt.contourf(lon, lat, cur_var[c_time, c_slice, :, :])
                    ax.imshow(cur_var[c_time, c_slice, :, :], extent=self.getExtent(lat, lon))
                    c_title = F'{c_var_name} {title} Z-level:{c_slice} Time:{c_time} '
                    plt.title(c_title, fontsize=20)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                plt.tight_layout()
                self._close_figure()

    def plot_4d_data_ncds(self, ds, var_names: list, z_levels: list, timesteps: list, title='',
                          file_name_prefix='', cmap='viridis'):
        """
        This is the main function to plot multiple z_levels from a NetCDF file (4D data)
        """
        checkFolder(self._output_folder)
        for c_slice in z_levels:
            for c_time in timesteps:
                plt.subplots(1, len(var_names), squeeze=True, figsize=(8*len(var_names), 8))
                for idx_var, c_var_name in enumerate(var_names):
                    cur_var = ds.variables.get(c_var_name)
                    ax = plt.subplot(1, len(var_names), idx_var+1)
                    plot_slice_eoa(cur_var[c_time, c_slice,:,:], ax, cmap=cmap)
                    c_title = F'{c_var_name} {title} Z-level:{c_slice} Time:{c_time} '
                    plt.title(c_title, fontsize=20)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()
    # ====================================== 3D Data =====================================

    def plot_3d_data_xarray_map(self, xr_ds, var_names: list, timesteps: list, title='', file_name_prefix='',
                                proj=ccrs.PlateCarree()):
        """
        Plots multiple z_levels from a NetCDF file (3D data). It plots the results in a map.
        It is assuming that the 3rd
        http://xarray.pydata.org/en/stable/plotting.html#maps
        """
        projection = proj

        checkFolder(self._output_folder)
        for c_time_idx in timesteps:
            plt.subplots(1, len(var_names), squeeze=True, figsize=(8*len(var_names), 8))
            for idx_var, c_var_name in enumerate(var_names):
                print(c_var_name)
                cur_var = xr_ds[c_var_name]
                # TODO. Hardcoded order Assuming the order of the coordinsates is lat, lon, time
                cur_coords_names = list(cur_var.coords.keys())
                # Assuming the order of the dims are time, lat, lon
                cur_dims_names = list(cur_var.dims)
                # TODO here the selection of the 'time' coordinate is hardcoded and changes depending
                # the file. Not sure how to fix it, so it will need to change for every problem. 
                # timevar = cur_coords_names[2]
                # Time selection 
                timevar = cur_coords_names[0] 
                c_time = xr_ds[timevar].values[c_time_idx]
                # Obtains data to the nearest requested time
                try:
                    cur_var = xr_ds[c_var_name].sel(**{timevar: c_time}, method='nearest')
                except Exception as e:
                    print(F"Warning for {c_var_name}!! (couldn't interpolate to the proper 'time' value: {e}")
                    cur_var = xr_ds[c_var_name].sel(**{cur_dims_names[0]: c_time_idx})
                ax = plt.subplot(1, len(var_names), idx_var+1, projection=projection)
                # ------------------ MAP STUFF -------------------------
                # https://rabernat.github.io/research_computing_2018/maps-with-cartopy.html
                cur_var.plot(ax=ax, transform=projection, cbar_kwargs={'shrink': 0.4})
                # ax.stock_img()  # Draws a basic topography of the world
                ax.coastlines(resolution='50m')  # Draws the coastline
                # ax = self.add_states(ax)
                ax = self.add_roads(ax)
                ax.add_feature(cartopy.feature.OCEAN)
                ax.add_feature(cartopy.feature.LAND, edgecolor='black')
                ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
                ax.add_feature(cartopy.feature.RIVERS)
                ax.gridlines()
                # ------------------ MAP STUFF -------------------------
                c_title = F'{c_var_name} {title} Time:{c_time} '
                plt.title(c_title, fontsize=20)

            file_name = F'{file_name_prefix}_{c_time}'
            pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
            self._close_figure()

    def plot_3d_data_ncds_map(self, ds, var_names: list, timesteps: list, title='',
                              file_name_prefix='', cmap='viridis', proj=ccrs.PlateCarree(), lonvar='Longitude',
                              latvar='Latitude'):
        """
        Plots multiple z_levels from a NetCDF file (4D data). It plots the results in a map.
        """
        lon = ds.variables[lonvar][:]
        lat = ds.variables[latvar][:]
        projection = proj

        checkFolder(self._output_folder)
        for c_time in timesteps:
            plt.subplots(1, len(var_names), squeeze=True, figsize=(8 * len(var_names), 8))
            for idx_var, c_var_name in enumerate(var_names):
                cur_var = ds.variables.get(c_var_name)
                ax = plt.subplot(1, len(var_names), idx_var + 1, projection=projection)
                # ------------------ MAP STUFF -------------------------
                # https://rabernat.github.io/research_computing_2018/maps-with-cartopy.html
                ax.stock_img()  # Draws a basic topography
                ax.coastlines(resolution='50m')  # Draws the coastline
                # ax = self.add_states(ax)
                ax = self.add_roads(ax)
                ax.add_feature(cartopy.feature.OCEAN)
                ax.add_feature(cartopy.feature.LAND, edgecolor='black')
                ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
                ax.add_feature(cartopy.feature.RIVERS)
                ax.gridlines()
                # ------------------ MAP STUFF -------------------------
                # plt.contourf(lon, lat, cur_var[c_time, c_slice, :, :])
                ax.imshow(cur_var[c_time, :, :], extent=self.getExtent(lat, lon))
                c_title = F'{c_var_name} {title} Time:{c_time} '
                plt.title(c_title, fontsize=20)

            file_name = F'{file_name_prefix}_{c_time:04d}'
            pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
            plt.tight_layout()

            self._close_figure()

    def plot_3d_data_ncds(self, ds, var_names: list, z_levels: list, title='',
                          file_name_prefix='', cmap='viridis'):
        """
        This is the main function to plot multiple z_levels (no time) from netCDF files
        """
        checkFolder(self._output_folder)
        for c_slice in z_levels:
                plt.subplots(1, len(var_names), squeeze=True, figsize=(8*len(var_names), 8))
                for idx_var, c_var_name in enumerate(var_names):
                    cur_var = ds.variables.get(c_var_name)
                    ax = plt.subplot(1, len(var_names), idx_var+1)
                    plot_slice_eoa(cur_var[c_slice,:,:], ax, cmap=cmap)
                    c_title = F'{c_var_name} {title} Z-level:{c_slice}'
                    plt.title(c_title, fontsize=20)
                plt.title("TEST", fontsize=30)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()


    def plot_3d_data_singlevar_np(self, data:list, z_levels: list, title='',
                          file_name_prefix='', cmap='viridis', max_imgs_per_row=4):
        """
        Plot 3D data from numpy arrays
        """
        checkFolder(self._output_folder)
        plt.subplots(1, len(z_levels), squeeze=True, figsize=(8 * max_imgs_per_row, 8*int(np.ceil(len(z_levels)/max_imgs_per_row))))
        for slice_idx, c_slice in enumerate(z_levels):
            # cur_img_row = np.floor(slice_idx/max_imgs_per_row) + 1
            # cur_img_col = slice_idx % max_imgs_per_row + 1
            ax = plt.subplot(1, len(z_levels), slice_idx+1)
            plot_slice_eoa(data[c_slice,:,:], ax, cmap=cmap)
            c_title = F'{title} Z-level:{c_slice}'
            plt.title(c_title, fontsize=20)

        file_name = F'{file_name_prefix}_{c_slice:04d}'
        # pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()


    def plot_3d_data_np(self, np_variables:list, var_names:list, z_levels: list, title='',
                          file_name_prefix='', cmap='viridis'):
        """
        Plot multiple z_levels.
        """
        checkFolder(self._output_folder)
        for c_slice in z_levels:
                plt.subplots(1, len(var_names), squeeze=True, figsize=(8*len(var_names), 8))
                for idx_var, c_var in enumerate(np_variables):
                    ax = plt.subplot(1, len(var_names), idx_var+1)
                    plot_slice_eoa(c_var[c_slice,:,:], ax, cmap=cmap)
                    if var_names != '':
                        c_title = F'{var_names[idx_var]} {title} Z-level:{c_slice}'
                    else:
                        c_title = F'{idx_var} {title} Z-level:{c_slice}'
                    plt.title(c_title, fontsize=20)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()


    def plot_1d_data_np(self, X, Ys,  title='', labels=[], file_name_prefix='', wide_ratio=1):
        """

        """
        plt.figure(figsize=[16*wide_ratio, 8])
        checkFolder(self._output_folder)
        for i, y in enumerate(Ys):
            if len(labels) > 0:
                assert len(labels) == len(Ys)
                plt.plot(X, y, self._COLORS[i], label=labels[i])
            else:
                plt.plot(X, y, self._COLORS[i])

        if len(labels) > 0:
            plt.legend()

        plt.title(title, fontsize=20)

        file_name = F'{file_name_prefix}'
        pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()
