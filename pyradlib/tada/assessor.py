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
import os
import pandas as pd
import pydicom
from . import complexity

# define default parameters
_default_integration_limit = 1
_default_integration_epsilon = 0.05
_default_identifies = ['name', 'mu']
_default_metrics = ['MCS']
_default_plan_prefix = 'rp'
_default_plan_extension = '.dcm'


def analyse_plan(path: str, metrics=_default_metrics):
    """ Calculate the complexity of a treatment plan.

    Args:
        path (str): The path of the radiotherapy treatment plan to be read.
        metrics (list): List of metrics to calculate.
    """
    ds = pydicom.read_file(path)
    columns, complexity_data = complexity.plan_complexity(ds)
    return pd.DataFrame(complexity_data, columns=columns)


def analyse_plan_archive(path: str, metrics=_default_metrics):
    """ Calculate the complexity for an archive of treatment plans.

    Args:
        path (str): The path of the radiotherapy treatment plan to be read.
        metrics (list): List of metrics to calculate.
    """
    complexity_data = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if name.lower().startswith(_default_plan_prefix) and name.lower().endswith(_default_plan_extension):
                ds = pydicom.read_file(os.path.join(root,name))
                columns, plan_complexity = complexity.plan_complexity(ds)
                complexity_data.extend(plan_complexity)
    return pd.DataFrame(complexity_data, columns = columns)