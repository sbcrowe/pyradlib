# -*- coding: utf-8 -*-
""" Radiochromic film calibration example module.

This module provides examples of application of film calibration functionality.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Emma Spelleken']
__license__ = 'GPL3'

# import required code
import radlib.film.calibration as cal

# define dose values
dose_extended = [0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 4000]
dose_standard = dose_extended[:10]

# define netod values (taken from Epson V800 and 10000XL scanners)
netod_extended_800 = [0.004751, 0.015599, 0.051076, 0.090765, 0.129687, 0.157541, 0.183196, 0.208629, 0.228482,
                      0.251974, 0.368132, 0.442065, 0.495993, 0.538633, 0.572647, 0.603522, 0.627520, 0.672461]
netod_standard_800 = netod_extended_800[:10]
netod_extended_10000 = [0.005468, 0.015911, 0.061723, 0.110433, 0.152401, 0.187004, 0.219768, 0.250146, 0.275071,
                        0.303447, 0.446947, 0.540865, 0.608908, 0.656808, 0.697727, 0.731365, 0.759067, 0.800121]
netod_standard_10000 = netod_extended_10000[:10]

# perform fits for standard dose range
standard_fits_800 = [cal.PolynomialFit(dose_standard, netod_standard_800, 4),
                     cal.NonLinearFit(dose_standard, netod_standard_800),
                     cal.DevicFit(dose_standard, netod_standard_800),
                     cal.TamponiFit(dose_standard, netod_standard_800)]
extended_fits_800 = [cal.PolynomialFit(dose_extended, netod_extended_800, 4),
                     cal.NonLinearFit(dose_extended, netod_extended_800),
                     cal.DevicFit(dose_extended, netod_extended_800),
                     cal.TamponiFit(dose_extended, netod_extended_800)]
standard_fits_10000 = [cal.PolynomialFit(dose_standard, netod_standard_10000, 4),
                       cal.NonLinearFit(dose_standard, netod_standard_10000),
                       cal.DevicFit(dose_standard, netod_standard_10000),
                       cal.TamponiFit(dose_standard, netod_standard_10000)]
extended_fits_10000 = [cal.PolynomialFit(dose_extended, netod_extended_10000, 4),
                       cal.NonLinearFit(dose_extended, netod_extended_10000),
                       cal.DevicFit(dose_extended, netod_extended_10000),
                       cal.TamponiFit(dose_extended, netod_extended_10000)]

# print data
print("Method, RMSE, CVRMSE, Mean relative error")
for fit in standard_fits_800:
    print(fit.name + ", " + str(fit.root_mean_square_error()) + " cGy, " +
          str(fit.coefficient_variation_root_mean_square_error() * 100) + "%, " +
          str(fit.mean_relative_error() * 100) + "%")
    fit.plot_calibration()
for fit in extended_fits_800:
    print(fit.name + ", " + str(fit.root_mean_square_error()) + " cGy, " +
          str(fit.coefficient_variation_root_mean_square_error() * 100) + "%, " +
          str(fit.mean_relative_error() * 100) + "%")
    fit.plot_calibration()
for fit in standard_fits_10000:
    print(fit.name + ", " + str(fit.root_mean_square_error()) + " cGy, " +
          str(fit.coefficient_variation_root_mean_square_error() * 100) + "%, " +
          str(fit.mean_relative_error() * 100) + "%")
    fit.plot_calibration()
for fit in extended_fits_10000:
    print(fit.name + ", " + str(fit.root_mean_square_error()) + " cGy, " +
          str(fit.coefficient_variation_root_mean_square_error() * 100) + "%, " +
          str(fit.mean_relative_error() * 100) + "%")
    fit.plot_calibration()