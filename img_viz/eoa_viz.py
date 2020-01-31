from os import listdir
from os.path import join

import matplotlib.pyplot as plt
from img_viz.common import *
from img_viz.constants import SliceMode, PlaneTypes, PlotMode
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

class EOAImageVisualizer:
    """This class makes plenty of plots from netcdf datasets or numpy arrays """
    _COLORS = ['y', 'r', 'c', 'b', 'g', 'w', 'k', 'y', 'r', 'c', 'b', 'g', 'w', 'k']
    _figsize = 8
    _font_size = 30
    _fig_prop = 1.8  # Proportion of each figure w/h
    _units = ''
    _max_imgs_per_row = 4

    def __init__(self, **kwargs):
        # All the arguments that are passed to the constructor of the class MUST have its name on it.
        self._disp_images = True
        self._output_folder = 'output'
        for arg_name, arg_value in kwargs.items():
            self.__dict__["_" + arg_name] = arg_value

    def __getattr__(self, attr):
        '''Generic getter for all the properties of the class'''
        return self.__dict__["_" + attr]

    def __setattr__(self, attr, value):
        '''Generic setter for all the properties of the class'''
        self.__dict__["_" + attr] = value

    def add_colorbar(self, fig, im, ax, show_color_bar):
        if show_color_bar:
            font_size_cbar = self._font_size * .6
            cbar = fig.colorbar(im, ax=ax)
            cbar.ax.tick_params(labelsize=font_size_cbar)
            cbar.set_label(self._units, fontsize=font_size_cbar*1.2)

    def get_proper_size(self, rows, cols):
        """
        Obtains the proper size for a figure.
        :param rows: how many rows will the figure have
        :param cols: how many colswill the figure have
        :param prop: Proportion is the proportion to use w/h
        :return:
        """
        if rows == 1:
            return self._figsize * cols * self._fig_prop, self._figsize
        else:
            return self._figsize * cols * self._fig_prop, self._figsize * rows

    def _close_figure(self):
        """Depending on what is disp_images, the figures are displayed or just closed"""
        if self._disp_images:
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
        create_folder(self._output_folder)
        for c_slice in z_levels:
            for c_time in timesteps:
                fig, axs = plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))
                for idx_var, c_var_name in enumerate(var_names):
                    # Depth selection, not sure if the name is always the same
                    c_depth = ds[z_levelvar].values[c_slice]
                    # Interpolates data to the nearest requested depth
                    cur_var = ds[c_var_name].sel(**{z_levelvar: c_depth}, method='nearest')
                    ax = plt.subplot(1, len(var_names), idx_var+1, projection=projection)
                    # ------------------ MAP STUFF -------------------------
                    # https://rabernat.github.io/research_computing_2018/maps-with-cartopy.html
                    cur_var.plot(ax=ax, transform=projection)
                    # ax.stock_img()  # Draws a basic topography of the world
                    ax.coastlines(resolution='50m')  # Draws the coastline
                    # ax = self.add_states(ax)
                    ax = self.add_roads(ax)
                    ax.add_feature(cartopy.feature.OCEAN)
                    ax.add_feature(cartopy.feature.LAND, edgecolor='black')
                    ax.add_feature(cartopy.feature.LAKES, edgecolor='black')
                    ax.add_feature(cartopy.feature.RIVERS)
                    # ------------------ MAP STUFF -------------------------
                    c_title = F'{c_var_name} {title} Z-level:{c_slice} Time:{c_time} '
                    ax.set_title(c_title, fontsize=self._font_size)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()

    def plot_4d_data_ncds_map(self, ds, var_names: list, z_levels: list, timesteps: list, title='',
                          file_name_prefix='', cmap='viridis', proj=ccrs.PlateCarree(),
                              lonvar='Longitude', latvar='Latitude', show_color_bar=True):
        """
        Plots multiple z_levels from a NetCDF file (4D data). It plots the results in a map.
        """
        lon = ds.variables[lonvar][:]
        lat = ds.variables[latvar][:]
        projection = proj

        create_folder(self._output_folder)
        for c_slice in z_levels:
            for c_time in timesteps:
                fig, axs = plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))
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
                    # ------------------ MAP STUFF -------------------------
                    # plt.contourf(lon, lat, cur_var[c_time, c_slice, :, :])
                    im = ax.imshow(cur_var[c_time, c_slice, :, :], extent=self.getExtent(lat, lon))
                    c_title = F'{c_var_name} {title} Z-level:{c_slice} Time:{c_time} '
                    ax.set_title(c_title, fontsize=self._font_size)
                    self.add_colorbar(fig, im, ax, show_color_bar)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                plt.tight_layout()
                self._close_figure()

    def plot_4d_data_ncds(self, ds, var_names: list, z_levels: list, timesteps: list, title='',
                          file_name_prefix='', cmap='viridis', show_color_bar=True):
        """
        This is the main function to plot multiple z_levels from a NetCDF file (4D data)
        """
        create_folder(self._output_folder)
        for c_slice in z_levels:
            for c_time in timesteps:
                fig, axs = plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))
                for idx_var, c_var_name in enumerate(var_names):
                    cur_var = ds.variables.get(c_var_name)
                    ax = plt.subplot(1, len(var_names), idx_var+1)
                    im = plot_slice_eoa(cur_var[c_time, c_slice,:,:], ax, cmap=cmap)
                    c_title = F'{c_var_name} {title} Z-level:{c_slice} Time:{c_time} '
                    ax.set_title(c_title, fontsize=self._font_size)
                    self.add_colorbar(fig, im, ax, show_color_bar)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()

    # ====================================== 3D Data =====================================

    def plot_3d_data_xarray_map(self, xr_ds, var_names: list, timesteps: list, title='', file_name_prefix='',
                                proj=ccrs.PlateCarree(), timevar_name='time'):
        """
        Plots multiple z_levels from a NetCDF file (3D data). It plots the results in a map.
        It is assuming that the 3rd
        http://xarray.pydata.org/en/stable/plotting.html#maps
        """
        projection = proj

        create_folder(self._output_folder)
        for c_time_idx in timesteps:
            plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))
            for idx_var, c_var_name in enumerate(var_names):
                print(c_var_name)
                cur_var = xr_ds[c_var_name]
                # TODO. Hardcoded order Assuming the order of the coordinsates is lat, lon, time
                cur_coords_names = list(cur_var.coords.keys())
                # Assuming the order of the dims are time, lat, lon
                cur_dims_names = list(cur_var.dims)
                c_time = xr_ds[timevar_name].values[c_time_idx]
                # Obtains data to the nearest requested time
                try:
                    cur_var = xr_ds[c_var_name].sel(**{timevar_name: c_time}, method='nearest')
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
                plt.title(c_title, fontsize=self._font_size)

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

        create_folder(self._output_folder)
        for c_time in timesteps:
            plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))
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
                plt.title(c_title, fontsize=self._font_size)

            file_name = F'{file_name_prefix}_{c_time:04d}'
            pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
            plt.tight_layout()

            self._close_figure()

    def plot_3d_data_ncds(self, ds, var_names: list, z_levels: list, title='',
                          file_name_prefix='', cmap='viridis'):
        """
        This is the main function to plot multiple z_levels (no time) from netCDF files
        """
        create_folder(self._output_folder)
        for c_slice in z_levels:
                plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))
                for idx_var, c_var_name in enumerate(var_names):
                    cur_var = ds.variables.get(c_var_name)
                    ax = plt.subplot(1, len(var_names), idx_var+1)
                    plot_slice_eoa(cur_var[c_slice,:,:], ax, cmap=cmap)
                    c_title = F'{c_var_name} {title} Z-level:{c_slice}'
                    plt.title(c_title, fontsize=self._font_size)
                plt.title("TEST", fontsize=30)

                file_name = F'{file_name_prefix}_{c_slice:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()


    def plot_3d_data_singlevar_np(self, data:list, z_levels= [], title='',
                          file_name_prefix='', cmap='viridis', flip_data=False,
                                  show_color_bar=True):
        """
        Plots all the z-layers for a single 3d var
        """
        create_folder(self._output_folder)
        rows = int(np.ceil(len(z_levels)/self._max_imgs_per_row))
        cols = int(min(self._max_imgs_per_row, len(z_levels)))
        if len(z_levels) == 0:
            z_levels = np.arange(data.shape[0])

        fig, axs = plt.subplots(rows, cols, squeeze=True,
                                figsize=self.get_proper_size(rows, cols))

        for slice_idx, c_slice in enumerate(z_levels):
            ax = plt.subplot(rows, cols, slice_idx+1)
            if flip_data:
                im = plot_slice_eoa(np.flip(np.flip(data[c_slice,:,:]),axis=1), ax, cmap=cmap)
            else:
                im = plot_slice_eoa(data[c_slice, :, :], ax, cmap=cmap)
            c_title = F'{title} Z-level:{c_slice}'
            ax.set_title(c_title, fontsize=self._font_size)
            self.add_colorbar(fig, im, ax, show_color_bar)

        file_name = F'{file_name_prefix}_{c_slice:04d}'
        pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()


    def plot_3d_data_np(self, np_variables:list, var_names:list, z_levels= [], title='',
                          file_name_prefix='', cmap='viridis', z_lavels_names = [], flip_data=False,
                        plot_mode=PlotMode.RASTER):
        """
        Plot multiple z_levels.
        """
        create_folder(self._output_folder)

        # If the user do not requires any z-leve, then all are plotted
        if len(z_levels) == 0:
            z_levels = range(np_variables[0].shape[0])

        for c_slice in z_levels:
                plt.subplots(1, len(var_names), squeeze=True, figsize=self.get_proper_size(1, len(var_names)))

                # Verify the index of the z_levels are the original ones.
                if len(z_lavels_names) != 0:
                    c_slice_txt = z_lavels_names[c_slice]
                else:
                    c_slice_txt = c_slice

                for idx_var, c_var in enumerate(np_variables):
                    ax = plt.subplot(1, len(var_names), idx_var+1)
                    if flip_data:
                        plot_slice_eoa(np.flip(np.flip(c_var[c_slice, :, :]),axis=1), ax, cmap=cmap, mode=plot_mode)
                    else:
                        plot_slice_eoa(c_var[c_slice,:,:], ax, cmap=cmap, mode=plot_mode)

                    if var_names != '':
                        c_title = F'{var_names[idx_var]} {title} Z-level:{c_slice_txt}'
                    else:
                        c_title = F'{idx_var} {title} Z-level:{c_slice_txt}'

                    plt.title(c_title, fontsize=self._font_size)

                file_name = F'{file_name_prefix}_{c_slice_txt:04d}'
                pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
                self._close_figure()

    # ====================================== 2D Data =====================================
    def plot_2d_data_np(self, np_variables: list, var_names: list, title='',
                        file_name_prefix='', cmap='viridis', flip_data=False,
                        plot_mode=PlotMode.RASTER, show_color_bar=True):
        """
        Plots an array of 2D fields that come as np arrays
        :param np_variables:
        :param var_names:
        :param title:
        :param file_name_prefix:
        :param cmap:
        :param flip_data:
        :param plot_mode:
        :return:
        """
        create_folder(self._output_folder)
        fig, axs = plt.subplots(squeeze=True, figsize=self.get_proper_size(1, len(var_names)), ncols=len(var_names))

        for idx_var, c_var in enumerate(np_variables):
            ax = plt.subplot(1, len(var_names), idx_var+1)
            if flip_data:
                im =plot_slice_eoa(np.flip(np.flip(c_var),axis=1), ax, cmap=cmap, mode=plot_mode)
            else:
                im =plot_slice_eoa(c_var, ax, cmap=cmap, mode=plot_mode)
            self.add_colorbar(fig, im, ax, show_color_bar)

            if var_names != '':
                c_title = F'{var_names[idx_var]} {title} '
            else:
                c_title = F'{idx_var} {title}'
            ax.set_title(c_title, fontsize=self._font_size)

        file_name = F'{file_name_prefix}'
        pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()

    # ====================================== 1D Data =====================================

    def plot_1d_data_np(self, X, Ys,  title='', labels=[], file_name_prefix='', wide_ratio=1):
        """

        """
        plt.figure(figsize=[16*wide_ratio, self._figsize])
        create_folder(self._output_folder)
        for i, y in enumerate(Ys):
            if len(labels) > 0:
                assert len(labels) == len(Ys)
                plt.plot(X, y, self._COLORS[i], label=labels[i])
            else:
                plt.plot(X, y, self._COLORS[i])

        if len(labels) > 0:
            plt.legend()

        plt.title(title, fontsize=self._font_size)

        file_name = F'{file_name_prefix}'
        pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()


    def plot_1d_data_xarray(self, xr_ds, var_names: list, title='', file_name_prefix=''):
        """
        Plots multiple variables from a NetCDF file (1D data). It plots the results in a line
        """
        create_folder(self._output_folder)
        for idx_var, c_var_name in enumerate(var_names):
            xr_ds[c_var_name].to_dataframe().plot()
            c_title = F'{c_var_name} {title}'
            plt.title(c_title, fontsize=self._font_size)

        file_name = F'{file_name_prefix}'
        pylab.savefig(join(self._output_folder, F'{file_name}.png'), bbox_inches='tight')
        self._close_figure()

