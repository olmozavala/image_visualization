from os import listdir
from os.path import join

from imageio import imwrite
import SimpleITK as sitk
from PIL import Image

from img_viz.common import *
from img_viz.constants import SliceMode, PlaneTypes


class MaskRCNNVisualizer:
    _disp_images = True
    _output_folder = 'output'
    _input_folder = 'input'
    _COLORS = ['y', 'r', 'c', 'b', 'g', 'w', 'k', 'y', 'r', 'c', 'b', 'g', 'w', 'k']

    _use_flipped_image = True

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

    def plot_for_mask_r_cnn_input(self, img_name, ctr_names,
                                  slices=SliceMode.ALL, plane=PlaneTypes.AXIAL,
                                  plot_empty_ctr=False):
        """
        Creates masks and 2d images from the original MRI. It stores them in the specific format for Mask R-CNN network.
        :param orig_slices:
        :param plot_empty_ctr:
        :param _use_flipped_image:
        :param input_img_folder:
        :param _output_folder:
        :param img_names:
        :param ctr_names:
        :return:
        """
        all_cases = listdir(self.input_folder)
        all_cases.sort()
        checkFolder(self.output_folder)

        for c_folder in all_cases:
            print(F'{join(self._output_folder, c_folder)}')

            # =========== Reading data ================
            temp = sitk.ReadImage(join(self.input_folder, c_folder, img_name))
            img_np = norm_image(sitk.GetArrayFromImage(temp))

            ctrs_np = []
            orig_ctr_names = ctr_names.copy()
            for id_ctr_name, c_ctr_name in enumerate(orig_ctr_names):  # Adding all ctrs (masks)
                try:
                    temp = sitk.ReadImage(join(self.input_folder, c_folder, c_ctr_name))
                    ctrs_np.append(sitk.GetArrayFromImage(temp))
                except Exception as e:
                    ctr_names.pop(id_ctr_name)
                    print(F'-----Warning for {c_ctr_name} (skiping this contour) error: {e} ----------------')

            c_slices = get_slices(slices, img_np, plane=plane)
            for c_slice in c_slices:
                # Indicates if we want to plot all the slices or only the one with contours
                if should_display_slice(ctrs_np, c_slice, plane, draw_only_ctrs=not(plot_empty_ctr)):

                    if self._use_flipped_image:
                        c_folder_name = F'{c_folder}_{plane.value}_{c_slice:04d}_f'
                        final_img_np = np.flip(get_proper_plane(img_np, plane, c_slice))
                    else:
                        c_folder_name = F'{c_folder}_{plane.value}_{c_slice:04d}'
                        final_img_np = get_proper_plane(img_np, plane, c_slice)

                    # Defining and creating folder for this image and slice
                    final_img_folder = join(self._output_folder, c_folder_name, 'images')
                    final_masks_folder = join(self._output_folder, c_folder_name, 'masks')
                    img_final_name = join(final_img_folder, F'{c_folder_name}.png')
                    checkFolder(final_img_folder)
                    checkFolder(final_masks_folder)

                    # Saving the imagme
                    imwrite(img_final_name, final_img_np)
                    print(F'---- {img_final_name} ---')

                    # ========= Saving each ctr as a binary mask ===========
                    for id_c_ctr, c_ctr_name in enumerate(ctr_names):
                        if self._use_flipped_image:
                            c_ctr_np = get_proper_plane(ctrs_np[id_c_ctr], plane, c_slice)
                        else:
                            c_ctr_np = np.flip(get_proper_plane(ctrs_np[id_c_ctr], plane, c_slice),1)
                        if np.sum(c_ctr_np) > 0:
                            # Save the ctr in the proper folder
                            ctr_final_name = join(final_masks_folder, F'{c_ctr_name}.png')
                            imwrite(ctr_final_name, 255 * c_ctr_np)
                            print(F'---- {ctr_final_name} ---')


    # def plot_multiple_histograms(all_hist, labels, save_file='', start_at=0, width=4):
    #     plt.figure(figsize=(12, 8))
    #     try:
    #         for ii, c_hist in enumerate(all_hist):
    #             x = c_hist[1][start_at:-1]
    #             y = c_hist[0][start_at:]
    #             plt.bar(x, y, width=width * np.ones(len(x)), alpha=.5, label=labels[ii])
    #
    #         plt.legend(loc='best')
    #
    #         if save_file != '':
    #             pylab.savefig(save_file, bbox_inches='tight')
    #
    #         _close_figure()
    #
    #     except Exception as e:
    #         print('---------------------------- Failed {} error: {} ----------------'.format(save_file, e))
    #
    # def plot_multiple_scatter(all_scatter, labels, save_file='', title=''):
    #     plt.figure(figsize=(12, 8))
    #     try:
    #         for ii, c_data in enumerate(all_scatter):
    #             c_data.sort()
    #             y = c_data
    #             plt.scatter(range(len(y)), y, label=labels[ii])
    #
    #         if title != '':
    #             plt.title(title)
    #         plt.legend(loc='best')
    #
    #         if save_file != '':
    #             pylab.savefig(save_file, bbox_inches='tight')
    #
    #         _close_figure()
    #     except Exception as e:
    #         print('---------------------------- Failed {} error: {} ----------------'.format(save_file, e))
    #
    # def plot_histograms_from_img(images_itk, title='', mode='2d', save_file=''):
    #     """
    #     Plots the histogram and one slice of N number of images
    #     :param save_file:
    #     :param mode:
    #     :param title: title of the figure
    #     :param images_itk:
    #     :return:
    #     """
    #     n_images = len(images_itk)
    #     fig, ax = plt.subplots(2, n_images, squeeze=True, figsize=(8 * n_images, 8))
    #     for ii, c_img in enumerate(images_itk):
    #         c_img = sitk.GetArrayFromImage(c_img)
    #         slices = get_slices(SliceMode.MIDDLE, np.array([c_img]))
    #
    #         if n_images > 1:
    #             t1 = ax[0][ii]
    #             t2 = ax[1][ii]
    #         else:
    #             t1 = ax[0]
    #             t2 = ax[1]
    #
    #         if mode == '2d':
    #             f_img = c_img[slices[0], :, :]
    #         else:
    #             f_img = c_img
    #
    #         t1.hist(f_img.flatten(), 'auto')
    #         t2.imshow(c_img[slices[0], :, :])
    #
    #     plt.title(title)
    #     if save_file != '':
    #         pylab.savefig(save_file, bbox_inches='tight')
    #
    #     _close_figure()



    # def draw_img_and_ctrs_3d_np(c_img: np.array, np_ctrs: list, slices=SliceMode.ALL, out_folder='', title='',
    #                             labels=[], draw_only_ctrs=True, all_dims=[0, 1, 2]) -> None:
    #     """
    #     This method draws a single image with multiple contours in a single figure.
    #     :param c_img: A 3D numpy array
    #     :param np_ctrs: A list of 3D numpy arrays with all the contour images.
    #     :param slices: One of the SliceModes or a numpy array with the slices indexes
    #     :param out_folder: Folder where to save the images.
    #     :param title: Used to display as title of the displayed image
    #     :param draw_only_ctrs: If true, only slices with contours will be stored/displayed.
    #     :param dim: the dimension to use to generate the plots  (-1 will draw images in every dimension)
    #     :param labels: list of strings with the labels used for each contour
    #     """
    #
    #     # Init figure
    #     tot_dims = len(all_dims)
    #
    #     # Iterate for each dimension
    #     for cur_dim_idx in range(tot_dims):
    #         cur_dim = all_dims[cur_dim_idx]
    #         cur_dim_size = c_img.shape[cur_dim_idx]
    #         slices = _get_slices(c_img, slices, cur_dim)
    #         # Iterate over each slice
    #         for cur_slice in slices:
    #             if _should_display_slice(np_ctrs, cur_slice, cur_dim, draw_only_ctrs):
    #                 try:
    #                     curax = plt.subplot(1, 1, 1)
    #                     curax.axis('off')
    #                     curax.imshow(c_img.take(cur_slice, cur_dim), cmap='gray')
    #                     # plt.colorbar(imres,ax=curax)
    #                     for idx_ctr, c_np_ctr in enumerate(np_ctrs):
    #                         cs = curax.contour(c_np_ctr.take(cur_slice, cur_dim), colors=COLORS[idx_ctr],
    #                                            linewidths=.4)
    #                         if len(labels) > 0:  # Adding the legends
    #                             curax.clabel(cs, inline=1, fontsize=0)
    #                             cs.collections[0].set_label(labels[idx_ctr])
    #
    #                     if len(labels) > 0:  # Position the legend
    #                         curax.legend(loc='upper right', framealpha=1, prop={'size': 15})
    #
    #                     c_title = F'{title} D{cur_dim}-S{cur_slice:03d}/{cur_dim_size} '
    #                     plt.title(c_title, fontsize=20)
    #                     file_name = F'D{cur_dim}_{cur_slice:03d}'
    #                     if out_folder != '':
    #                         if not Path(out_folder).exists():
    #                             os.makedirs(str(out_folder))
    #                         pylab.savefig(join(out_folder, F'{file_name}.jpg'), bbox_inches='tight')
    #                     _close_figure()
    #
    #                 except Exception as e:
    #                     print(
    #                         F'---------------ERROR Failed for D{cur_dim}-S{cur_slice:03d}/{cur_dim_size} error: {e} ----------------')

