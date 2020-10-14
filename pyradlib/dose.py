# -*- coding: utf-8 -*-
""" Dose module.

This module provides a class for encapsulation of n-dimension dose distributions.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__license__ = "GPL3"

# import required code
import numpy as np
import scipy.ndimage as ndimage
from bisect import bisect


class Dose:
    """Encapsulates dose distribution data."""

    def __init__(self, dose_data, coordinate_data):
        """Instantiates the calibration fit using a non-linear function.

        Args:
            dose_data (np.array): The calibration doses delivered.
            coordinate_data (np.array): The coordinates of the element boundaries.
        """
        # assert pre-conditions
        assert len(dose_data.shape) == len(coordinate_data), 'Dose and coordinate dimensionality do not match.'
        # instantiate class variables
        self.dose_data = dose_data
        self.coordinate_data = coordinate_data

    def dose(self, coordinate):
        """Determines dose at nominal coordinate, by interpolation or nearest neighbour.

        Args:
            coordinate (tuple): The coordinate of the dose to be determined (in mm).

        Returns:
            dose (float): The dose (in Gy) at the specified coordinate.

        Raises:
            ValueError: If number of coordinate points incorrect.
        """
        # pre-conditions
        if not len(coordinate) == len(self.coordinate_data):
            raise ValueError('Number of coordinate points incorrect.')
        # determine index and return dose
        index_coord = self.index_coordinate(self, coordinate)
        return ndimage.map_coordinates(self.dose_data, index_coord, mode='nearest')

    def index_coordinate(self, coordinate):
        """Determines the (float) index coordinates corresponding to a physical coordinate.

        Args:
            coordinate (tuple): The physical coordinate of the point to be converted (in mm).

        Returns:
            The index coordinate corresponding to the physical coordinate.

        Raises:
            ValueError: If number of coordinate points incorrect.
        """
        # pre-conditions
        if not len(coordinate) == len(self.coordinate_data):
            raise ValueError('Number of coordinate points incorrect.')
        # determine indices
        index_coordinates = []
        for dimension in len(coordinate):
            dimension_coords = self.coordinate_data[dimension]
            if not self.increasing_coordinates(dimension):
                dimension_coords = list(reversed(dimension_coords))
            base_index = bisect(dimension_coords, coordinate[dimension])
            offset = 0
            # assess extrapolation cases
            if coordinate[dimension] < dimension_coords[0]:
                offset = 0
            elif coordinate[dimension] >= dimension_coords[-1]:
                offset = 0
            else:
                offset = ((coordinate[dimension] - self.coordinate_data[base_index])
                          / (self.coordinate_data[base_index + 1] - self.coordinate_data[base_index]))
            dimension_index = base_index + offset
            if not self.increasing_coordinates(dimension):
                dimension_index = len(dimension_coords) - dimension_index
            index_coordinates.add(dimension_index)
        return tuple(index_coordinates)

    def increasing_coordinates(self, dimension):
        """Determines if coordinates in dimension increase with increasing indices.

        Args:
            dimension (int): The dimension to be evaluated (0-based).

        Returns:
            bool: True if coordinates in dimension are increasing, False otherwise.

        Raises:
            IndexError: If specified dimension dose not exist.
            ValueError: If specified dimension has only one discrete coordinate.
        """
        # pre-conditions
        if not 0 <= dimension < len(self.coordinate_data):
            raise IndexError('Specified dimension does not exist.')
        if not len(self.coordinate_data[dimension] > 1):
            raise ValueError('Specified dimension has only one discrete coordinate.')
        # return result
        return self.coordinate_data[dimension][0] < self.coordinate_data[dimension][1]