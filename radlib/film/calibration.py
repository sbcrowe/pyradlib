# -*- coding: utf-8 -*-
""" Radiochromic film calibration module.

This module provides dose calibration functionality for radiochromic film.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Tanya Kairn', 'Samuel Peet']
__license__ = "GPL3"

# import required code
import numpy as np
from scipy import optimize
from numpy.polynomial.polynomial import polyval

# define default parameters
_default_number_format = '{0:0.2f}'


class NonLinearFit:
    """Performs calibration fit using Devic non-linear fitting function."""

    name = 'Non-linear fit'
    variable = 'netOD'

    def __init__(self, dose_data, netod_data):
        """Instantiates the calibration fit using a non-linear function.

        Args:
            dose_data (np.array): The calibration doses delivered.
            netod_data (np.array): The net optical densities of irradiated film.
        """
        # assert pre-conditions
        assert len(dose_data) == len(netod_data), 'Number of dose and netOD data points disagree.'
        assert len(dose_data) > 4, 'Insufficient number of data points for non-linear fit.'
        # instantiate class variables
        self.dose_data = dose_data
        self.netod_data = netod_data
        self.minimum = np.min(netod_data)
        self.maximum = np.max(netod_data)
        popt, pcov = optimize.curve_fit(self._non_linear_fitting_function,
                                        netod_data, dose_data, p0=[1, 0, 2])
        self.a, self.b, self.n = popt
        perr = np.sqrt(np.diag(pcov))
        self.a_error, self.b_error, self.n_error = perr
        self.degrees_of_freedom = len(dose_data) - 4
        calculated_dose = self.dose(np.array(netod_data))
        self.rmse = np.sqrt(np.sum((np.array(dose_data) - calculated_dose) ** 2) / self.degrees_of_freedom)
        self.cvrmse = self.rmse / np.mean(dose_data)

    def _non_linear_fitting_function(x, a, b, n):
        """The Devic non-linear fitting function."""
        return a * x ** n + b * x

    def dose(self, netod):
        """Returns dose calculated for provided netOD."""
        return self._non_linear_fitting_function(netod, self.a, self.b, self.n)

    def error(self, netod, netod_error):
        """Returns dose error for provided netOD."""
        A = netod**self.n * self.a_error
        B = netod * self.b_error
        C = self.a * netod**self.n * np.log(self.n) * self.n_error
        D = (self.a * self.n * netod**(self.n-1) + self.b) * netod_error
        return np.sqrt(A**2 + B**2 + C**2 + D**2)

    def description(self):
        """Returns text describing the calibration fit."""
        return ('d(x)=' + _default_number_format.format(self.a) + 'x^'
                + _default_number_format.format(self.n) + ' + '
                + _default_number_format.format(self.b) + 'x + '
                + _default_number_format.format(self.c))


class PolynomialFit:
    """Performs calibration fit using polynomial fitting function."""

    name = 'Polynomial fit'
    variable = 'netOD'

    def __init__(self, dose_data, netod_data, degree):
        """Instantiates the calibration fit using a polynomial.
        
        Args:
            dose_data (np.array): The calibration doses delivered.
            netod_data (np.array): The net optical densities of irradiated film.
            degree: The degree of the fitting polynomial.
        """
        # assert pre-conditions
        assert len(dose_data) == len(netod_data), 'Number of dose and netOD data points disagree.'
        assert len(dose_data) > (degree + 1), 'Insufficient number of data points for polynomial fit.'
        # instantiate class variables
        self.dose_data = dose_data
        self.netod_data = netod_data
        self.minimum = np.min(netod_data)
        self.maximum = np.max(netod_data)
        self.degree = degree
        self.coefficients, covariance_matrix = np.polyfit(netod_data, dose_data, degree, cov=True)
        self.errors = np.sqrt(np.diag(covariance_matrix))
        self.degrees_of_freedom = len(dose_data) - (degree + 1)
        calculated_dose = self.dose(np.array(netod_data))
        self.rmse = np.sqrt(np.sum((np.array(dose_data) - calculated_dose) ** 2) / self.degrees_of_freedom)
        self.cvrmse = self.rmse / np.mean(dose_data)

    def dose(self, netod):
        """Returns dose corresponding to provided netOD."""
        return polyval(netod, self.coefficients)

    def error(self, netod, netod_error):
        """Returns dose error for provided netOD."""
        sum_square_fitting_parameters = 0
        multiplier_for_noise = 0
        for degree_index in np.arange(0, self.degree):
            power = self.degree - degree_index
            coeff = self.coefficients[degree_index]
            error = self.errors[degree_index]
            sum_square_fitting_parameters += (netod ** power * error) ** 2
            if (power > 0):
                multiplier_for_noise += coeff * power * netod**(power - 1)
        variance_noise = (multiplier_for_noise * netod_error)**2
        return np.sqrt(sum_square_fitting_parameters + variance_noise)

    def description(self):
        """Returns text describing the calibration fit."""
        text = 'd(x)='
        for index in np.arange(len(self.coef)):
            text += (' + ' + _default_number_format.format(self.coef[1])
                     + 'x^' + str(index+1))
        return text.replace('x^1','x').replace('x^0','').replace('= +','=')


class TamponiFit:
    """Performs calibration fit using Tamponi fitting function."""

    name = 'Tamponi fit'
    variable = 'netOD'

    def __init__(self, dose_data, netod_data):
        """Instantiates the calibration fit using the Tamponi method."""
        # assert pre-conditions
        assert len(dose_data) == len(netod_data), 'Number of dose and netOD data points disagree.'
        assert len(dose_data) > 1, 'Insufficient number of data points for Tamponi fit.'
        # instantiate class variables
        self.dose_data = dose_data
        self.netod_data = netod_data
        self.minimum = np.min(netod_data)
        self.maximum = np.max(netod_data)
        self.ref_dose = dose_data[-1]
        self.ref_netod = netod_data[-1]
        R = np.array(dose_data) / self.ref_dose
        eR = np.array(netod_data) / self.ref_netod
        popt, pcov = optimize.curve_fit(self._tamponi_fitting_function, eR, R, p0=[0.5])
        self.a = popt[0]
        calculated_dose = self.dose(np.array(netod_data))
        self.degrees_of_freedom = len(dose_data) - 1
        self.rmse = np.sqrt(np.sum((np.array(dose_data) - calculated_dose) ** 2) / self.degrees_of_freedom)
        self.cvrmse = self.rmse / np.mean(dose_data)

    def _tamponi_fitting_function(R, a):
        """The Tamponi fitting function."""
        return np.log10(1 + (a - 1) * R) / np.log10(a)

    def dose(self, netod):
        """Returns dose corresponding to provided netOD."""
        eR = netod / self.ref_netod
        # return self.ref_dose * np.log(1+(self.a-1)*eR) / np.log(self.a)
        return self._tamponi_fitting_function(eR, self.a)

    def description(self):
        """Returns text describing the calibration fit."""
        return 'a=' + _default_number_format.format(self.a)


class LewisFit:
    """Performs calibration fit using Lewis fitting function."""

    name = 'Lewis fit'
    variable = 'Net response'

    def __init__(self, dose_data, response_data):
        """Instantiates the calibration fit using the Tamponi method."""
        # assert pre-conditions
        assert len(dose_data) == len(response_data), 'Number of dose and net response data points disagree.'
        assert len(dose_data) > 3, 'Insufficient number of data points for Lewis fit.'
        # instantiate class variables
        self.dose_data = dose_data
        self.response_data = response_data
        self.minimum = np.min(response_data)
        self.maximum = np.max(response_data)
        popt, pcov = optimize.curve_fit(self._lewis_fitting_function, response_data, dose_data, p0=[0, 1, 0])
        self.a, self.b, self.c = popt
        calculated_dose = self.dose(np.array(response_data))
        self.degrees_of_freedom = len(dose_data) - 3
        self.rmse = np.sqrt(np.sum((np.array(dose_data) - calculated_dose) ** 2) / self.degrees_of_freedom)
        self.cvrmse = self.rmse / np.mean(dose_data)

    def _lewis_fitting_function(x, a, b, c):
        """The Lewis fitting function."""
        return ((b) / (x + a)) + c

    def dose(self, response):
        """Returns dose corresponding to provided net response."""
        return self._lewis_fitting_function(response, self.a, self.b, self.c)

    def description(self):
        """Returns text describing the calibration fit."""
        return ('d(x)=((' + _default_number_format.format(self.b)
                + ') / (x+' + _default_number_format.format(self.a)
                + ')) + ' + _default_number_format.format(self.c))


class YaoFit:
    """Performs calibration fit using Yao fitting function.

    Note: the equation presented by Yao: PR = (1+hD/m) / (1+D/m) can be
    inverted as D = m(PR-1) / (h-PR); where PR = PV / PV_0
    """

    name = 'Yao fit'
    variable = 'Net response'

    def __init__(self, dose_data, response_data):
        """Instantiates the calibration fit using the Yao method."""
        # assert pre-conditions
        assert len(dose_data) == len(response_data), 'Number of dose and net response data points disagree.'
        assert len(dose_data) > 2, 'Insufficient number of data points for Yao fit.'
        # instantiate class variables
        self.dose_data = dose_data
        self.response_data = response_data
        self.minimum = np.min(response_data)
        self.maximum = np.max(response_data)
        median_dose = np.median(dose_data)
        popt, pcov = optimize.curve_fit(self._yao_fitting_function, response_data, dose_data, p0=[0.19, median_dose])
        self.h, self.m = popt
        calculated_dose = self.dose(np.array(response_data))
        self.degrees_of_freedom = len(dose_data) - 2
        self.rmse = np.sqrt(np.sum((np.array(dose_data) - calculated_dose) ** 2) / self.degrees_of_freedom)
        self.cvrmse = self.rmse / np.mean(dose_data)

    def _yao_fitting_function(x, a, b):
        """The Yao fitting function."""
        return (b * (x - 1) / (a - x))

    def dose(self, response):
        """Returns dose corresponding to provided response ratio."""
        return self._yao_fitting_function(response, self.h, self.m)

    def description(self):
        """Returns text describing the calibration fit."""
        return ('d(x)= ' + _default_number_format.format(self.m)
                + '(x - 1) / (' + _default_number_format.format(self.h)
                + ' - x)')