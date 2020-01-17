from enum import Enum


class SliceMode(Enum):
    """ Enum to select the different types of slices"""
    ALL = 1
    MIDDLE = 2
    MIDDLE_THIRD = 3


class PlaneTypes(Enum):
    AXIAL = 'AX'
    CORONAL = 'COR'
    SAGITTAL = 'SAG'
    ALL = 'ALL'

class PlotMode(Enum):
    """ Enum to select the different types of slices"""
    RASTER=1
    CONTOUR=2
    MERGED=3

