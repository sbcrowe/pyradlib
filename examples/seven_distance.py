# -*- coding: utf-8 -*-
""" Seven-distance calibration example module.

This module provides examples of application of the seven distance method for brachytherapy source calibration.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Emily Simpson-Page']
__license__ = 'GPL3'

# add parent directory to path, in case it is not already
import sys
import os
from pathlib import Path
sys.path.append(os.fspath(Path(__file__).parent.parent.resolve()))

# import required code
import pyradlib.brachy.calibration as cal
from rich.console import Console
from rich.table import Table
console = Console()

# instantiate measurement (nC), distance (cm) pair array
_data = [[2.34239e-8, 0.16], [1.51045e-8, 0.2], [1.05189e-8, 0.24], [7.7443e-9, 0.28], [5.98325e-9, 0.32], [4.75618e-9, 0.36], [3.88959e-9, 0.4]]
_methods = {'Least squares': cal.whole_leastsquares, 'Piecewise fsolve': cal.piecewise_fsolve}

table = Table(show_header=True)
table.add_column("Method")
table.add_column("F")
table.add_column("Ms")
table.add_column("c")

# perform calculations
for method in _methods:
    f, ms, c = _methods[method](_data)
    table.add_row(method, str(f), str(ms), str(c))
    f, ms, c = _methods[method](_data, once_removed=True)
    table.add_row(method + '*', str(f), str(ms), str(c))
console.print(table)