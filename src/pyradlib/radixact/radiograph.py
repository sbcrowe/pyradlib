# -*- coding: utf-8 -*-
"""Radiograph analysis module.

This module provides functionality for processing of radiograph data from Synchrony treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import os

import numpy as np
import numpy.typing as npt


class RadixactRadiograph:
    # region Constructors

    def __init__(self, image: npt.ArrayLike, uid: str = None) -> RadixactRadiograph:
        """Initialises the radiograph image class.

        Parameters
        ----------
        image : npt.ArrayLike
            Numpy array containing pixel values.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

        Returns
        -------
        RadixactRadiograph
            The radiograph image, encapsulated in a helper class.
        """
        self._image = image
        self.uid = uid

    @classmethod
    def from_bin(cls, path: str | os.PathLike, uid: str = None) -> RadixactRadiograph:
        """Reads radiograph image from a binary file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the binary file.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

        Returns
        -------
        RadixactRadiograph
            The radiograph image, encapsulated in a helper class.

        Note
        ----
        Binary files containing radiographs or digitally reconstructed radiographs
        consist of 960x960 16-bit signed integers (big-endian). For each radiograph
        a pixel represents 9 clustered detector elements (i.e., 3x3 elements from the
        2880x2880 element array).
        """
        data = np.fromfile(path, dtype=np.int16).reshape((960, 960))
        return cls(data, uid)

    # endregion
