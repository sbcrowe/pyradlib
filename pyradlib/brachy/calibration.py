# -*- coding: utf-8 -*-
""" Brachytherapy air-kerma strength calibration module.

This module provides functions assisting in the determination of brachytherapy source air-kerma strength.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Emily Simpson-Page']
__license__ = 'GPL3'

# import required code
import numpy as np
import itertools 
from scipy.optimize import fsolve, least_squares
from sympy.core.symbol import symbols
from sympy.solvers.solveset import nonlinsolve
from gekko import GEKKO

# define global variables
_max_c = 0.01

def piecewise_fsolve(data, average='mean', start_fmc=None, once_removed=False):
    """ Solves overdetermined system of non-linear equations by piecewise averaging of fsolve solutions.

    Args:
        data (array): array of arrays containing measurement and distance pairs
        average (str): Mode for averaging, 'mean' or 'median'.
        start_fmc (array): array of starting optimisation values for f, Ms and c.
        once_removed (bool): flag for calculation of f using averaged Ms and c values, instead of using the average f value.

    Returns:
        Array containing F, M_s and c results.
    """
    results = []
    min_md = np.min(np.array(data)[:,0])
    if start_fmc is None:
        start_fmc = np.array([min_md,min_md/2,0])
    for combination in itertools.combinations(data, 3):
        def func(z):
            f, m, c = z
            F = np.empty((3))
            F[0] = (combination[0][0]-m)*((c+combination[0][1])**2) - f
            F[1] = (combination[1][0]-m)*((c+combination[1][1])**2) - f
            F[2] = (combination[2][0]-m)*((c+combination[2][1])**2) - f
            return F
        result = fsolve(func,start_fmc)
        # sanity tests (primary + scatter > scatter, positive scatter)
        if (result[0] >= result[1]) and (result[1] >= 0) and (abs(result[2]) <= _max_c): 
            results.append(fsolve(func,start_fmc))
    if average is 'mean':
        return_result = np.average(np.array(results), axis=0)
    elif average is 'median':
        return_result = np.median(np.array(results), axis=0)
    if once_removed:
        f_values = (np.array(data)[:,0]-return_result[1])*(np.array(data)[:,1]+return_result[2])**2
        if average is 'mean':
            return_result[0] = np.mean(f_values)
        elif average is 'median':
            return_result[0] = np.median(f_values)
    return return_result

def piecewise_nonlinsolve(data, average='mean', once_removed=False):
    """ Solves overdetermined system of non-linear equations by piecewise averaging of nonlinsolve solutions.

    Args:
        data (array): array of arrays containing measurement and distance pairs
        average (str): Mode for averaging, 'mean' or 'median'.
        once_removed (bool): flag for calculation of f using averaged Ms and c values, instead of using the average f value.

    Returns:
        Array containing F, M_s and c results.
    """
    # note: this approach appears to be consistent with WolframAlpha approach.
    results = []
    for combination in itertools.combinations(data, 3):
        f, m, c = symbols('f, m, c', real=True)
        eq1 = (combination[0][0] - m) * ((combination[0][1] + c)**2) - f
        eq2 = (combination[1][0] - m) * ((combination[1][1] + c)**2) - f
        eq3 = (combination[2][0] - m) * ((combination[2][1] + c)**2) - f
        # catch exceptions and discard solutions that are not sensible.
        try: 
            sols = nonlinsolve([eq1,eq2,eq3], [f, m, c])
            for solution in sols:
                if (solution[0] >= solution[1]) and (solution[1] >=0) and (abs(solution[2]) <= _max_c):
                    results.append(solution)
        except:
            sols = []
    if average is 'mean':
        return_result = np.average(np.array(results), axis=0)
    elif average is 'median':
        return_result = np.median(np.array(results), axis=0)
    if once_removed:
        f_values = (np.array(data)[:,0]-return_result[1])*(np.array(data)[:,1]+return_result[2])**2
        if average is 'mean':
            return_result[0] = np.mean(f_values)
        elif average is 'median':
            return_result[0] = np.median(f_values)
    return return_result

def piecewise_gekko(data, average='mean', once_removed=False):
    """ Solves overdetermined system of non-linear equations by piecewise averaging of gekko solutions.

    Args:
        data (array): array of arrays containing measurement and distance pairs
        average (str): Mode for averaging, 'mean' or 'median'.
        once_removed (bool): flag for calculation of f using averaged Ms and c values, instead of using the average f value.

    Returns:
        Array containing F, M_s and c results.
    """
    results = []
    max_ms = np.min(np.array(data)[:,0])
    for combination in itertools.combinations(data, 3):
        model = GEKKO()
        f = model.Var(max_ms, 1e-12)
        m = model.Var(max_ms/2, 1e-12, max_ms)
        c = model.Var(0, -_max_c, _max_c)
        model.Equation((combination[0][0]-m)*(c+combination[0][1])**2-f==0)
        model.Equation((combination[1][0]-m)*(c+combination[1][1])**2-f==0)
        model.Equation((combination[2][0]-m)*(c+combination[2][1])**2-f==0)
        model.solve(disp=False)
        # discard solutions that are not sensible
        solution = [f.value[0], m.value[0], c.value[0]]
        if (solution[0] >= solution[1]) and (solution[1] >=0) and (abs(solution[2]) <= _max_c) and (1e-12 not in solution):
            results.append(solution)
    if average is 'mean':
        return_result = np.average(np.array(results), axis=0)
    elif average is 'median':
        return_result = np.median(np.array(results), axis=0)
    if once_removed:
        f_values = (np.array(data)[:,0]-return_result[1])*(np.array(data)[:,1]+return_result[2])**2
        if average is 'mean':
            return_result[0] = np.mean(f_values)
        elif average is 'median':
            return_result[0] = np.median(f_values)
    return return_result

def whole_leastsquares(data, start_fmc=None, once_removed=False):
    """ Solves overdetermined system of non-linear equations by least squares fitting.

    Args: 
        data (array): array of arrays containing measurement and distance pairs
        start_fmc (array): array of initial F, M_s and c estimates.
        once_removed (bool): flag for calculation of f using averaged Ms and c values, instead of using the average f value

    Returns:
        Array containing F, M_s and c results.
    """
    min_md = np.min(np.array(data)[:,0])
    if start_fmc is None:
        start_fmc = [min_md, min_md/3, 0]
    bounds = ([0,0,-0.01],[min_md, min_md, 0.01])
    def func(args):
        r = (np.array(data)[:,0]-args[1])*((args[2]+np.array(data)[:,1])**2) - args[0]
        return r
    optimise_result = least_squares(func, start_fmc, jac='3-point', bounds=bounds, method='trf', ftol=None, xtol=1e-12, gtol=None, tr_solver='lsmr')
    if (once_removed):
        return np.array([np.mean((np.array(data)[:,0]-optimise_result.x[1])*(np.array(data)[:,1]-optimise_result.x[2])**2), optimise_result.x[1], optimise_result.x[2]])
    else:
        return optimise_result.x
