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

    def __init__(self, image: npt.ArrayLike) -> RadixactRadiograph:
        """Initialises the radiograph image class.

        Parameters
        ----------
        image : npt.ArrayLike
            Numpy array containing pixel values.

        Returns
        -------
        RadixactRadiograph
            The radiograph image, encapsulated in a helper class.
        """
        self._image = image

    @classmethod
    def from_bin(cls, path: str | os.PathLike) -> RadixactRadiograph:
        """Reads radiograph image from a binary file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the binary file.

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
        return cls(data)

    @classmethod
    def from_compressed(cls, path: str | os.PathLike) -> RadixactRadiograph:
        """Reads radiograph image from an npz file.

        Parameters
        ----------
        path : str | os.PathLake
            Path to the npz file.

        Returns
        -------
        RadixactRadiograph
            The radiograph image, encapsulated in a helper class.
        """
        cls.from_npz(path)

    @classmethod
    def from_npz(cls, path: str | os.PathLake) -> RadixactRadiograph:
        """Reads radiograph image from an npz file.

        Parameters
        ----------
        path : str | os.PathLake
            Path to the npz file.

        Returns
        -------
        RadixactRadiograph
            The radiograph image, encapsulated in a helper class.
        """
        with np.load(path) as npz_data:
            image = dict(npz_data)["image"]
            return cls(image)

    # endregion

    # region Public methods

    def to_compressed(self, path: str | os.PathLike) -> None:
        """Writes radiograph to compressed npz file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for npz file to be written.
        """
        self.to_npz(path)

    def to_npz(self, path: str | os.PathLike, compress: bool = True) -> None:
        """Writes motion data to numpy npz file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for npz file to be written.
        compress : bool, optional
            Flag to indicate whether to compress npz. Default is True.
        """
        data_dict = {}
        data_dict["image"] = self._image
        method = np.savez_compressed if compress else np.savez
        method(path, **data_dict)

    # endregion
