from os import listdir
from os.path import join

import matplotlib.pyplot as plt
import pylab
import scipy.misc as misc
import SimpleITK as sitk

from img_viz.common import *
from img_viz.constants import SliceMode, PlaneTypes


class MedicalImageVisualizer:
    _disp_images = True
    _output_folder = 'output'
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

    def _plot_imgs_and_ctrs_np(self, config) -> None:
        """
        This is the main function to plot multiple slices.
        """
        imgs = config['imgs']
        ctrs = config['ctrs']
        orig_slices = config['slices']
        plane = config['plane']
        title = config['title']
        draw_only_ctrs = config['draw_only_ctrs']
        file_name_prefix = config['file_name_prefix']
        labels = config['labels']

        checkFolder(self._output_folder)
        if plane != PlaneTypes.ALL:  # Draw single plane
            # This should plot a single image with all the contours overlayed
            first_img_shape = imgs[0].shape
            # Validate the size of the images
            for c_img in imgs:
                if c_img.shape != first_img_shape:
                    raise Exception('The shape of the images must be the same')

            slices = get_slices(orig_slices, imgs[0], plane)
            for c_slice in slices:
                draw_slice = should_display_slice(ctrs, c_slice, plane, draw_only_ctrs)
                if draw_slice:
                    fig, ax = plt.subplots(1, len(imgs), squeeze=True, figsize=(8, 8))
                    for img_idx, c_img in enumerate(imgs):
                        img_slice = get_proper_plane(c_img, plane, c_slice)
                        ctrs_slice = [get_proper_plane(np_ctr, plane, c_slice) for np_ctr in ctrs]
                        if len(imgs) > 1:
                            plot_slice(img_slice, ctrs_slice, ax[img_idx], labels)
                        else:
                            plot_slice(img_slice, ctrs_slice, ax, labels)
                        c_title = F'{title} {plane.value} {c_slice:04d}'
                        file_name = F'{file_name_prefix}_{plane.value}_{c_slice:04d}'
                        plt.title(c_title, fontsize=20)
                    pylab.savefig(join(self._output_folder, F'{file_name}.jpg'), bbox_inches='tight')
        else:
            if len(imgs) != 1:
                raise Exception('The number of image allowed for Plane type ALL must be 1')
            # In this case it should plot 3 images, one for each plane. Here we force it to plot
            # only the middle slice
            c_img = imgs[0]
            plt.subplots(1, 3, squeeze=True, figsize=(8 * 3, 8))
            for id_plane, plane in enumerate([PlaneTypes.AXIAL, PlaneTypes.SAGITTAL, PlaneTypes.CORONAL]):
                ax = plt.subplot(1, 3, id_plane + 1)
                c_slice = get_slices(SliceMode.MIDDLE, c_img, plane)[0]
                img_slice = get_proper_plane(c_img, plane, c_slice)
                ctrs_slice = [get_proper_plane(np_ctr, plane, c_slice) for np_ctr in ctrs]
                plot_slice(img_slice, ctrs_slice, ax, labels)
                c_title = F'{title} ALL {c_slice:04d}'
                plt.title(c_title, fontsize=20)

            file_name = F'{file_name_prefix}_{plane.value}_{c_slice:04d}'
            pylab.savefig(join(self._output_folder, F'{file_name}.jpg'), bbox_inches='tight')

        self._close_figure()

    def plot_imgs_and_ctrs_itk(self, itk_imgs: list, itk_ctrs: list, slices=SliceMode.ALL, title='',
                              plane=PlaneTypes.AXIAL, draw_only_ctrs=True, file_name_prefix='', labels=[]):
        """
        Main function to display images that come from an itk format
        :param itk_imgs : itk 3D image
        :param itk_ctrs: list of itk 3D contours
        :param slices: range of slices or SliceMode
        :param title: string title of the plot
        :param plane: PlaneTypes  Which plane to save
        :param draw_only_ctrs: boolean Indicates if we want to plot only slices where there is a contour
        :param file_name_prefix: A name to add in front of the file name
        :param labels: a list of strings to label each contour
        :return:
        """
        args = {
            'imgs': [sitk.GetArrayFromImage(itk_img) for itk_img in itk_imgs],
            'ctrs': [sitk.GetArrayFromImage(itk_ctr) for itk_ctr in itk_ctrs],
            'slices': slices,
            'plane': plane,
            'draw_only_ctrs': draw_only_ctrs,
            'file_name_prefix': file_name_prefix,
            'labels': labels,
            'title': title
        }
        self._plot_imgs_and_ctrs_np(args)

    def plot_img_and_ctrs_itk(self, itk_img, itk_ctrs: list, slices=SliceMode.ALL, title='',
                              plane=PlaneTypes.AXIAL, draw_only_ctrs=True, file_name_prefix='', labels=[]):
        """
        Main function to display images that come from an itk format
        :param itk_img : itk 3D image
        :param itk_ctrs: list of itk 3D contours
        :param slices: range of slices or SliceMode
        :param title: string title of the plot
        :param plane: PlaneTypes  Which plane to save
        :param draw_only_ctrs: boolean Indicates if we want to plot only slices where there is a contour
        :param file_name_prefix: A name to add in front of the file name
        :param labels: a list of strings to label each contour
        :return:
        """
        args = {
            'imgs': [sitk.GetArrayFromImage(itk_img)],
            'ctrs': [sitk.GetArrayFromImage(itk_ctr) for itk_ctr in itk_ctrs],
            'slices': slices,
            'plane': plane,
            'draw_only_ctrs': draw_only_ctrs,
            'file_name_prefix': file_name_prefix,
            'labels': labels,
            'title': title
        }
        self._plot_imgs_and_ctrs_np(args)

    def plot_locations_np(self, c_img, locations, file_name_prefix=''):
        """ This function plots the locations received in locations. For each plane.
        :param c_img:
        :param locations:
        :param file_name_prefix:
        """
        other_axes = [[1, 2], [2, 0], [1, 0]]
        shape_str = str(c_img.shape)
        checkFolder(self._output_folder)
        for id_loc, c_loc in enumerate(locations):
            plt.subplots(1, 3, squeeze=True, figsize=(8 * 3, 8))
            for id_plane, plane in enumerate([PlaneTypes.AXIAL, PlaneTypes.SAGITTAL, PlaneTypes.CORONAL]):
                c_loc = np.array(c_loc)
                c_dim = get_axis_idx(plane)
                c_slice = c_loc[c_dim]

                if c_slice < c_img.shape[c_dim]:
                    pt = c_loc[other_axes[c_dim]]
                    ax = plt.subplot(1, 3, id_plane + 1)
                    ax.imshow(get_proper_plane(c_img, plane, c_slice), cmap='gray')
                    ax.scatter(pt[0],pt[1], color=self._COLORS[id_loc])
                    c_title = F'{c_slice} {pt} -- {plane.value} -- {shape_str} '
                    plt.title(c_title, fontsize=20)
                else:
                    raise Exception(
                        F'The slice {c_slice} is larger than the shape of the array {c_img.shape[c_dim]}')

            pylab.savefig(join(self._output_folder, F'{file_name_prefix}_Loc_{id_loc}.jpg'), bbox_inches='tight')
            self._close_figure()

    def plot_multiple_bar_plots(self, tuple_dicts, file_name='', title='', legends=[]):
        """
        Plots the DSC as a bar plot
        :param tuple_dicts: List with multiple dictionaries with string keys (Case-0000) and values (DSC)
        :param file_name:
        :param legends:
        :param title:
        :return:
        """
        tot_columns = len(tuple_dicts)  # Number of Dice coefficients to compare
        tot_ex = len(tuple_dicts[0])  # Number of examples
        plt.figure(figsize=(8 * tot_ex * tot_columns / 14, 8))

        min_val = np.finfo('f').max
        max_val = np.finfo('f').min

        for cur_data in tuple_dicts:
            cur_min = min(cur_data.values())
            cur_max = max(cur_data.values())
            if cur_min < min_val:
                min_val = cur_min
            if cur_max > max_val:
                max_val = cur_max

        # Iterate over each column
        bar_width = .9
        for ii, cur_data in enumerate(tuple_dicts):
            # The trick here is to decide where to position each column
            x_axis = np.arange(ii, (tot_columns * tot_ex + ii), tot_columns)
            # print(F'Positions of each bar {x_axis}')
            if legends != '':
                plt.bar(x_axis, height=tuple_dicts[ii].values(),
                        width=bar_width,
                        tick_label=list(tuple_dicts[ii].keys()), align='edge', label=legends[ii])
            else:
                plt.bar(x_axis, tuple_dicts[ii].values(),
                        tick_label=list(tuple_dicts[ii].keys()), align='edge')
        # Defines where to put the X labels, position and labels
        x_axis_ticks  = np.arange( 0, tot_columns * tot_ex, tot_columns)
        plt.xticks(x_axis_ticks, list(tuple_dicts[0].keys()), rotation=20)
        # Defines the limits of the X and Y axis
        plt.ylim([min_val - .1, max(max_val + .1, 1)])
        plt.xlim([-1, tot_columns * tot_ex * 1.01])
        plt.legend(loc='best')
        plt.grid()
        plt.title(title)
        if file_name != '':
            plt.savefig(join(self._output_folder,file_name), bbox_inches='tight')
        self._close_figure()

    def plot_histograms_from_imgs(self, images_itk: list, file_name='', titles=[],
                                  plane=PlaneTypes.AXIAL):
        """
        This functino creates an image with two rows, in one row is an example middle slice
        of the image and in the other row the corresponding histogram. If PlaneType = ALL then the histogram
        is from the whole 3D image.
        :param plane:
        :param file_name:
        :param title: title of the figure
        :param images_itk:
        :return:
        """

        n_images = len(images_itk)
        fig, ax = plt.subplots(2, n_images, squeeze=True, figsize=(8 * n_images, 8))
        # Iterate over all the images
        for ii, c_img in enumerate(images_itk):
            img_np = sitk.GetArrayFromImage(c_img)
            middle_slice = get_slices(SliceMode.MIDDLE, img_np)[0]

            if n_images > 1:
                t1 = ax[0][ii]
                t2 = ax[1][ii]
            else:
                t1 = ax[0]
                t2 = ax[1]

            if plane == PlaneTypes.ALL:
                data_for_histogram = img_np
            else:
                data_for_histogram = get_proper_plane(img_np, plane, middle_slice)

            t1.hist(data_for_histogram.flatten(), 'auto')
            t2.imshow(get_proper_plane(img_np, plane, middle_slice))
            if len(titles) > 0:
                t1.title.set_text(titles[ii])
                t2.title.set_text(titles[ii])

        if file_name != '':
            pylab.savefig(join(self._output_folder, file_name), bbox_inches='tight')

        self._close_figure()
