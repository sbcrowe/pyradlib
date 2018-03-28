# -*- coding: utf-8 -*-
""" Radiochromic film calibration module.

This module provides dose calibration functionality for radiochromic film.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Tanya Kairn', 'Samuel Peet']
__license__ = 'GPL3'

# import required code
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

# define default parameters
_default_number_format = '{0:0.2f}'
_default_plot_resolution = 1000
_default_fill_alpha = 0.1
_default_max_fev = 10000


def plot_calibration(calibration, filename):
    """ Plot the calibration as .png file.
    
    Args:
        calibration: The calibration to be plotted.
        filename: The path of the file to be saved.
    """
    os.makedirs(filename, exist_ok=True)
    fitting_data_x = np.linspace(calibration.minimum, calibration.maximum, _default_plot_resolution)
    plt.plot(fitting_data_x, calibration.dose(fitting_data_x), 'b-')
    plt.fill_between(fitting_data_x, calibration.dose(fitting_data_x) + calibration.error(fitting_data_x),
                     calibration.dose(fitting_data_x) - calibration.error(fitting_data_x),
                     color='blue', alpha=_default_fill_alpha)
    plt.plot(calibration.netod_data, calibration.dose_data, 'b.')
    plt.xlabel(calibration.variable)
    plt.ylabel('Dose (Gy)')
    plt.title(calibration.name + '; ' + calibration.equation() + '; RMSE = ' + calibration.rmse)
    plt.savefig(filename)
    plt.close()


class Fit:
    """Parent class for calibration fitting functions."""

    def root_mean_square_error(self):
        """Determine the root mean square error of the calibration fit (in dose units)."""
        return np.sqrt(np.sum((self.dose_data - self.calculated_dose) ** 2) / self.degrees_of_freedom)

    def coefficient_variation_root_mean_square_error(self):
        """Determine the coefficient of variation of the root mean square error."""
        return self.root_mean_square_error() / np.mean(self.dose_data)

    def mean_relative_error(self):
        """Returns mean of relative errors between fitting model and known doses."""
        if 0 in self.dose_data:
            return np.sum(abs(self.dose_data[1:] - self.calculated_dose[1:]) / self.dose_data[1:]) \
                   / (len(self.dose_data)-1)
        else:
            return np.sum(abs(self.dose_data - self.calculated_dose) / self.dose_data) / len(self.dose_data)

    def akaike_information_criterion(self, effective_uncertainty):
        """Returns Akaike information criteria."""
        return (2 * self.fit_parameters) + np.sum((np.array(self.dose_data) - self.calculated_dose) ** 2) / (
            self.degrees_of_freedom * effective_uncertainty ** 2) + (
            2 * self.fit_parameters * (self.fit_parameters + 1) / (self.degrees_of_freedom - 1))

    def plot_calibration(self):
        film_data_range = np.linspace(np.min(self.film_data), np.max(self.film_data), _default_plot_resolution)
        plt.plot(film_data_range, self.dose(film_data_range), 'b-')
        plt.fill_between(film_data_range, self.dose(film_data_range) + self.root_mean_square_error(),
                         self.dose(film_data_range) - self.root_mean_square_error(), color='blue', alpha='0.5')
        plt.plot(self.film_data, self.dose_data, 'b.')
        plt.xlabel(self.variable)
        plt.ylabel('Dose (Gy)')
        plt.title(self.name + '; ' + self.description())
        plt.show()

class NonLinearFit(Fit):
    """Performs calibration fit using non-linear fitting function."""

    name = 'Non-linear fit'
    variable = 'netOD'
    fit_parameters = 3

    def __init__(self, dose_data, netod_data):
        """Instantiates the calibration fit using a non-linear function.

        Args:
            dose_data (np.array): The calibration doses delivered.
            netod_data (np.array): The net optical densities of irradiated film.
            
        Raises:
            ValueError: If number of dose and netOD data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(netod_data):
            raise ValueError('Number of dose and netOD data points disagree.')
        if not len(dose_data) > 4:
            raise ValueError('Insufficient number of data points for non-linear fit.')
        # instantiate class variables
        self.dose_data = np.array(dose_data)
        self.film_data = np.array(netod_data)
        popt, pcov = optimize.curve_fit(self._non_linear_fitting_function, netod_data, dose_data,
                                        p0=[5000, 1000, 2.5], maxfev=_default_max_fev)
        self.a, self.b, self.n = popt
        perr = np.sqrt(np.diag(pcov))
        self.a_error, self.b_error, self.n_error = perr
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters
        self.calculated_dose = self.dose(self.film_data)

    @staticmethod
    def _non_linear_fitting_function(x, a, b, n):
        """The Devic non-linear fitting function."""
        return a * x ** n + b * x

    def dose(self, netod):
        """Returns dose calculated for provided netOD."""
        return self._non_linear_fitting_function(netod, self.a, self.b, self.n)

    def error(self, netod, netod_error):
        """Returns dose error for provided netOD."""
        A = netod ** self.n * self.a_error
        B = netod * self.b_error
        C = self.a * netod ** self.n * np.log(self.n) * self.n_error
        D = (self.a * self.n * netod ** (self.n - 1) + self.b) * netod_error
        return np.sqrt(A ** 2 + B ** 2 + C ** 2 + D ** 2)

    def description(self):
        """Returns text describing the calibration fit."""
        return ('d(x)=' + _default_number_format.format(self.a) + 'x^'
                + _default_number_format.format(self.n) + ' + '
                + _default_number_format.format(self.b) + 'x')


class DevicFit(NonLinearFit):
    """Performs calibration fit using Devic non-linear fitting function."""

    name = 'Devic fit'
    fit_parameters = 2

    def __init__(self, dose_data, netod_data):
        """Instantiates the calibration fit using a Devic non-linear function.

        Args:
            dose_data (np.array): The calibration doses delivered.
            netod_data (np.array): The net optical densities of irradiated film.

        Raises:
            ValueError: If number of dose and netOD data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(netod_data):
            raise ValueError('Number of dose and netOD data points disagree.')
        if not len(dose_data) > 4:
            raise ValueError('Insufficient number of data points for non-linear fit.')
        # instantiate class variables
        self.dose_data = np.array(dose_data)
        self.film_data = np.array(netod_data)
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters
        # perform regression
        possible_n = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        best_rmse = float('inf')
        best_parameters = 0, 0, 0
        best_pcov = 0, 0
        for test_n in possible_n:
            test_popt, test_pcov = optimize.curve_fit(
                lambda x, a, b: self._non_linear_fitting_function(x, a, b, test_n), netod_data, dose_data,
                p0=[5000, 1000], maxfev=_default_max_fev)
            test_a, test_b = test_popt
            test_dose = self._non_linear_fitting_function(self.film_data, test_a, test_b, test_n)
            test_rmse = np.sqrt(np.sum((self.dose_data - test_dose) ** 2) / self.degrees_of_freedom)
            if test_rmse < best_rmse:
                best_rmse = test_rmse
                best_parameters = test_a, test_b, test_n
                best_pcov = test_pcov
        self.a, self.b, self.n = best_parameters
        perr = np.sqrt(np.diag(best_pcov))
        self.a_error, self.b_error = perr
        self.calculated_dose = self.dose(self.film_data)


class PolynomialFit(Fit):
    """Performs calibration fit using polynomial fitting function."""

    name = 'Polynomial fit'
    variable = 'netOD'

    def __init__(self, dose_data, netod_data, degree):
        """Instantiates the calibration fit using a polynomial.
        
        Args:
            dose_data (np.array): The calibration doses delivered.
            netod_data (np.array): The net optical densities of irradiated film.
            degree: The degree of the fitting polynomial.
            
        Raises:
            ValueError: If number of dose and netOD data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(netod_data):
            raise ValueError('Number of dose and netOD data points disagree.')
        if not len(dose_data) > degree + 1:
            raise ValueError('Insufficient number of data points for polynomial fit.')
        # instantiate class variables
        self.dose_data = np.array(dose_data)
        self.film_data = np.array(netod_data)
        self.degree = degree
        self.coefficients, covariance_matrix = np.polyfit(netod_data, dose_data, degree, cov=True)
        self.errors = np.sqrt(np.diag(covariance_matrix))
        self.fit_parameters = self.degree + 1
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters
        self.calculated_dose = self.dose(np.array(netod_data))

    def dose(self, netod):
        """Returns dose corresponding to provided netOD."""
        return np.polyval(self.coefficients, netod)

    def error(self, netod, netod_error):
        """Returns dose error for provided netOD."""
        sum_square_fitting_parameters = 0
        multiplier_for_noise = 0
        for degree_index in np.arange(0, self.degree):
            power = self.degree - degree_index
            coeff = self.coefficients[degree_index]
            error = self.errors[degree_index]
            sum_square_fitting_parameters += (netod ** power * error) ** 2
            if power > 0:
                multiplier_for_noise += coeff * power * netod ** (power - 1)
        variance_noise = (multiplier_for_noise * netod_error) ** 2
        return np.sqrt(sum_square_fitting_parameters + variance_noise)

    def description(self):
        """Returns text describing the calibration fit."""
        text = 'd(x)='
        for index in np.arange(len(self.coefficients)):
            text += (' + ' + _default_number_format.format(self.coefficients[index])
                     + 'x^' + str(index + 1))
        return text.replace('x^1', 'x').replace('x^0', '').replace('= +', '=')


class TamponiFit(Fit):
    """Performs calibration fit using Tamponi fitting function."""

    name = 'Tamponi fit'
    variable = 'netOD'
    fit_parameters = 1

    def __init__(self, dose_data, netod_data):
        """Instantiates the calibration fit using the Tamponi method.

        Args:
            dose_data (np.array): The calibration doses delivered.
            netod_data (np.array): The net optical densities of irradiated film.

        Raises:
            ValueError: If number of dose and netOD data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(netod_data):
            raise ValueError('Number of dose and netOD data points disagree.')
        if not len(dose_data) > 1:
            raise ValueError('Insufficient number of data points for Tamponi fit.')
        # instantiate class variables
        self.dose_data = np.array(dose_data)
        self.film_data = np.array(netod_data)
        self.ref_dose = self.dose_data[-1]
        self.ref_netod = self.film_data[-1]
        R = np.array(dose_data) / self.ref_dose
        eR = np.array(netod_data) / self.ref_netod
        popt, pcov = optimize.curve_fit(self._tamponi_fitting_function, eR, R, p0=[0.5])
        self.a = popt[0]
        self.calculated_dose = self.dose(np.array(netod_data))
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters

    @staticmethod
    def _tamponi_fitting_function(R, a):
        """The Tamponi fitting function."""
        return np.log10(1 + (a - 1) * R) / np.log10(a)

    def dose(self, netod):
        """Returns dose corresponding to provided netOD."""
        eR = netod / self.ref_netod
        return self.ref_dose * self._tamponi_fitting_function(eR, self.a)

    def description(self):
        """Returns text describing the calibration fit."""
        return 'a=' + _default_number_format.format(self.a)


class LewisFit(Fit):
    """Performs calibration fit using Lewis fitting function."""

    name = 'Lewis fit'
    variable = 'Net response'
    fit_parameters = 3

    def __init__(self, dose_data, response_data):
        """Instantiates the calibration fit using the Lewis method.
        
        Args:
            dose_data (np.array): The calibration doses delivered.
            response_data (np.array): The net response values of irradiated film.
        
        Raises:
            ValueError: If number of dose and net response data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(response_data):
            raise ValueError('Number of dose and net response data points disagree.')
        if not len(dose_data) > 3:
            raise ValueError('Insufficient number of data points for Lewis fit.')
        # instantiate class variables
        self.dose_data = np.array(dose_data)
        self.film_data = np.array(response_data)
        popt, pcov = optimize.curve_fit(self._lewis_fitting_function, response_data, dose_data, p0=[0, 1, 0])
        self.a, self.b, self.c = popt
        self.calculated_dose = self.dose(np.array(response_data))
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters

    @staticmethod
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


class YaoFit(Fit):
    """Performs calibration fit using Yao fitting function.

    Note: the equation presented by Yao: PR = (1+hD/m) / (1+D/m) can be
    inverted as D = m(PR-1) / (h-PR); where PR = PV / PV_0
    """

    name = 'Yao fit'
    variable = 'Net response'
    fit_parameters = 2

    def __init__(self, dose_data, response_data):
        """Instantiates the calibration fit using the Yao method.
        
        Args:
            dose_data (np.array): The calibration doses delivered.
            response_data (np.array): The net response values of irradiated film.
        
        Raises:
            ValueError: If number of dose and net response data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(response_data):
            raise ValueError('Number of dose and net response data points disagree.')
        if not len(dose_data) > 2:
            raise ValueError('Insufficient number of data points for Yao fit.')
        # instantiate class variables
        self.dose_data = np.array(dose_data)
        self.film_data = np.array(response_data)
        median_dose = np.median(dose_data)
        popt, pcov = optimize.curve_fit(self._yao_fitting_function, response_data, dose_data, p0=[0.19, median_dose])
        self.h, self.m = popt
        self.calculated_dose = self.dose(np.array(response_data))
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters

    @staticmethod
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


class SumSignalFit(Fit):
    """Performs calibration fit using sum signal fitting function, as suggested by ."""

    name = 'Sum signal fit'
    variable = 'Sum signal'
    fit_parameters = 4

    def __init__(self, dose_data, sum_signal_data):
        """Instantiates the calibration fit using the Yao method.

        Args:
            dose_data (np.array): The calibration doses delivered.
            sum_signal_data (np.array): The sum signal values of irradiated film.

        Raises:
            ValueError: If number of dose and sum signal data points incorrect.
        """
        # pre-conditions
        if not len(dose_data) == len(sum_signal_data):
            raise ValueError('Number of dose and sum signal data points disagree.')
        if not len(dose_data) > 4:
            raise ValueError('Insufficient number of data points for Sum Signal fit.')
        # instantiate class variables
        self.dose_data = dose_data
        self.film_data = sum_signal_data
        popt, pcov = optimize.curve_fit(self._sum_signal_fitting_function, sum_signal_data, dose_data, p0=[1, 1, 1, 1])
        self.a, self.b, self.c, self.d = popt
        self.calculated_dose = self.dose(np.array(sum_signal_data))
        self.degrees_of_freedom = len(dose_data) - self.fit_parameters

    @staticmethod
    def _sum_signal_fitting_function(x, a, b, c, d):
        """The sum signal double exponential fitting function."""
        return a * np.exp(b * x) + c * np.exp(d * x)

    def dose(self, sum_signal):
        """Returns dose corresponding to provided sum signal."""
        return self._sum_signal_fitting_function(sum_signal, self.a, self.b, self.c, self.d)

    def description(self):
        """Returns text describing the calibration fit."""
        return ('d(x)= ' + _default_number_format.format(self.a)
                + ' exp(' + _default_number_format.format(self.b)
                + 'x) + ' + _default_number_format.format(self.c)
                + ' exp(' + _default_number_format.format(self.d)
                + 'x)')