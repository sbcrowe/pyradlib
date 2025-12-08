# -*- coding: utf-8 -*-
""" Statistical process control module.

This module provides statistical process control methods to assess QA results.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__license__ = "GPL3"

# import required code
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# find data
_data_root_path = os.path.join(os.fspath(Path(__file__).parent.resolve()),'data')
_constants = pd.read_csv(os.path.join(_data_root_path,'constants.csv'), delimiter=',')

def tolerance_limit(x: list[float], m : int = 2, k : int = 3):
    """ Calculate tolerance limit based on data.

    Args:
        x (list[float]): The data to be analysed.
        m (int): The sample size for moving ranges.
        k (int): The number of sigma used as a limit.
    """
    constants = _constants[_constants['m'] == m]
    mR_x = np.abs(np.array(x[1:]) - np.array(x[:-1]))
    cl_x = np.mean(x)
    cl_r = np.mean(mR_x)
    sd_x = np.std(x)
    lcl_x = np.mean(x) - (np.mean(mR_x) * k / constants['d2'])


def xmr_chart(x: list[float], k : int = 3):
    """ Calculate tolerance limit based on data.

    Args:
        x (list[float]): The data to be analysed.
        k (int): The number of sigma used as a limit.
    """
    constants = _constants[_constants['m'] == 2]
    r = np.abs(np.array(x[1:]) - np.array(x[:-1]))
    cl_x = np.mean(x)
    cl_r = np.mean(r)
    sd_x = np.std(x)
    ucl_x = float(np.mean(x) + (np.mean(r) * k / constants['d2']))
    lcl_x = float(np.mean(x) - (np.mean(r) * k / constants['d2']))
    ucl_r = float(np.mean(r) * constants['D4'])
    lcl_r = float(np.mean(r) * constants['D3'])
    # Plot x and mR charts
    fig, axs = plt.subplots(2, figsize=(15,15), sharex=True)
    axs[0].plot(x, linestyle='-', marker='o', color='black')
    axs[0].axhline(cl_x, color='blue')
    axs[0].axhline(ucl_x, color = 'red', linestyle = 'dashed')
    axs[0].axhline(lcl_x, color = 'red', linestyle = 'dashed')
    axs[0].set_title('Individual Chart')
    axs[0].set(xlabel='Case', ylabel='Value')
    axs[1].plot(r, linestyle='-', marker='o', color='black')
    axs[1].axhline(cl_r, color='blue')
    axs[1].axhline(ucl_r, color='red', linestyle ='dashed')
    axs[1].axhline(lcl_r, color='red', linestyle ='dashed')
    axs[1].set_ylim(bottom=0)
    axs[1].set_title('mR Chart')
    axs[1].set(xlabel='Case', ylabel='Range')