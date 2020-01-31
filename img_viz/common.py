
import matplotlib.pyplot as plt
import pylab
from img_viz.constants import SliceMode, PlaneTypes, PlotMode
import numpy as np
import os

_COLORS = ['y', 'r', 'c', 'b', 'g', 'w', 'k', 'y', 'r', 'c', 'b', 'g', 'w', 'k']


def create_folder(output_folder):
    """ It only creates a folder if it doesn't exist"""
    if not(os.path.exists(output_folder)):
        os.makedirs(output_folder)

def plot_slice_eoa(c_img, ax, cmap='gray', mode=PlotMode.RASTER) -> None:
    """
    Plots an simple 2D img for EOA data.
    :param c_img: 2D Image
    :param ax:
    :return:
    """
    c_ax = ax
    c_ax.axis('off')

    if mode == PlotMode.RASTER:
        im = c_ax.imshow(c_img, cmap=cmap)
    if mode == PlotMode.CONTOUR:
        im = c_ax.contour(c_img)
    if mode == PlotMode.MERGED:
        im = c_ax.imshow(c_img, cmap=cmap)
        c_ax.contour(c_img, colors='r')

    return im


def plot_slice(c_img, ctrs: list, ax, labels, cmap='gray') -> None:
    """
    Plots an image with its contours.
    :param c_img: 2D Image
    :param ctrs:  Array of
    :param ax:
    :param labels: list of string with labels
    :return:
    """
    ctr_line_width = 1  # Larger is larger
    label_font_size = 15
    c_ax = ax
    c_ax.axis('off')
    c_ax.imshow(c_img, cmap=cmap)
    if len(ctrs) > 0:
        for idx, c_ctr in enumerate(ctrs):
            cs = c_ax.contour(c_ctr, colors=_COLORS[idx], linewidths=ctr_line_width)
            if len(labels) > 0:
                c_ax.clabel(cs, inline=1, fontsize=0)  # fsize = 0 avoids showing the contour numbers
                cs.collections[0].set_label(labels[idx])
        if len(labels) > 0:
            c_ax.legend(loc='upper right', framealpha=1, prop={'size': label_font_size})


def get_slices(orig_slices, arr, plane=PlaneTypes.AXIAL):
    """Get the slices for this array (normally is the number we specify), but
    we can ask for the middle slice or inside third of the slices
    :param PlaneTypes plane:
    :param SliceMode orig_slices: It can also be a numpy array
    :param arr: 3D numpy array
    """
    # Used to decide which axis to use to get the number of slices
    dim_idx = get_axis_idx(plane)
    cur_dim = arr.shape[dim_idx]
    # Deciding if we draw all the images or not
    if isinstance(orig_slices, SliceMode):
        if orig_slices == SliceMode.ALL:
            slices = range(cur_dim)
        elif orig_slices == SliceMode.MIDDLE:
            slices = [int(cur_dim / 2)]
        elif orig_slices == SliceMode.MIDDLE_THIRD:
            bottom = int(np.floor(cur_dim / 3))
            top = int(np.ceil(cur_dim * (2 / 3)))
            slices = range(bottom, top)
        else:
            raise Exception(F'The "slices" option is incorrect: {orig_slices}')
    else:
        slices = orig_slices
    return slices


def should_display_slice(ctrs, slice_num, plane, draw_only_ctrs=False):
    """Function that decide if we should display or not a slice"""
    # The idea is that we should display the slice if there is at least one contour inside it
    if (len(ctrs) == 0) or (not (draw_only_ctrs)):
        return True
    else:
        # Only draw slices where there is at least one contour
        if len(ctrs) > 0:
            for np_ctr in ctrs:
                c_ctr = get_proper_plane(np_ctr, plane, slice_num)
                if np.amax(c_ctr) > 0:  # Avoid black slices
                    return True

    return False  # In this case we didn't find any contour


def get_axis_idx(plane):
    """
    It indicates which index correspond to each plane
    """
    dim_idx = 0
    if plane == PlaneTypes.SAGITTAL:
        dim_idx = 2
    if plane == PlaneTypes.CORONAL:
        dim_idx = 1
    if plane == PlaneTypes.AXIAL:
        dim_idx = 0
    return dim_idx


def norm_image(img_np):
    """ Normalizes an image in the range 0 - 255 """
    return (img_np - np.amin(img_np))/np.ptp(img_np)


def get_proper_plane(arr, plane, slice_num):
    """
    This function selects the proper 2D slice from the slice number and plane axis
    :param arr: 3D numpy array
    :param plane: Enum with the plane
    :param slice: slice number
    :return:
    """
    cur_ax = get_axis_idx(plane)
    # Avoid index out of bounds for images with different # of slices
    if slice_num < arr.shape[cur_ax]:
        if cur_ax == 0:
            img_slice = arr[slice_num,:,:]
        if cur_ax == 1:
            img_slice = arr[:,slice_num,:]
        if cur_ax == 2:
            img_slice = arr[:, :, slice_num]
        return img_slice
    else:
        raise Exception(F'The slice {slice_num} is larger than the shape of the array {arr.shape[get_axis_idx(plane)]}')
