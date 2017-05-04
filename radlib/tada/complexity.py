# -*- coding: utf-8 -*-
""" Treatment complexity assessment module..

This module provides radiotherapy treatment complexity assessment functionality.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Tanya Kairn']
__license__ = "GPL3"

# import required code
import dicom
import numpy as np
from PIL import Image

# define default parameters
_default_fluence_bixel_number = 160
_default_integration_limit = 1
_default_integration_epsilon = 0.05


def save_fluence_image(fluence, filename):
    """ Write the fluence map as an 8-bit grayscale .png file.

    Args:
        fluence: The fluence map to be plotted.
        filename: The path of the file to be saved.
    """
    offset = -fluence.min()
    scale = 255 / (fluence.max() - fluence.min())
    grayscale_img = Image.fromarray(((fluence + offset) * scale).astype(np.uint8))
    grayscale_img.save(filename)
