from os import listdir
from os.path import join

import matplotlib.pyplot as plt
import pylab
import scipy.misc as misc
import SimpleITK as sitk

from img_viz.common import *
from img_viz.constants import SliceMode, PlaneTypes


class TimeSeriesVisualizer:
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

    def plot_multicolumns_from_df(self, df, column_names, x_axis='def', title='', legends=[],
                                  file_name='') -> None:
        """
        This is the main function to plot multiple slices.
        """
        create_folder(self._output_folder)

        if x_axis == 'def':
            x = np.arange(len(df.loc[:,column_names[0]].values))
        else:
            x = df.loc[:, x_axis].values

        fig, ax = plt.subplots()
        for i, cur_col in enumerate(column_names):
            cur_col_data = df.loc[:,cur_col].values

            if len(legends) > 0:
                ax.scatter(x, cur_col_data, c=self._COLORS[i], label=legends[i])
            else:
                ax.scatter(x, cur_col_data, c=self._COLORS[i], label=column_names[i])

        ax.legend()
        ax.grid(True)
        plt.xlabel('Observations')
        plt.ylabel('Model')
        plt.title(title, fontsize=20)
        min_val = min(x)
        max_val = max(x)
        plt.plot((min_val,max_val), (min_val, max_val), 'k')
        pylab.savefig(join(self._output_folder, F'{file_name}'), bbox_inches='tight')
        plt.show()

        self._close_figure()


