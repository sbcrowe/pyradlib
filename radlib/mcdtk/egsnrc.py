# -*- coding: utf-8 -*-
""" EGSnrc module.

This module provides EGSnrc file handling functionality.
"""

# define default parameters
_default_electron_energy_cutoff = 0.521
_default_photon_energy_cutoff = 0.01
_default_exact_boundary_crossing = True
_default_low_energy_simulation = False


def transport_parameter_text(electron_energy_cutoff=_default_electron_energy_cutoff,
                             photon_energy_cutoff=_default_photon_energy_cutoff,
                             exact_boundary_crossing=_default_exact_boundary_crossing,
                             low_energy_simulation=_default_low_energy_simulation):
    """ Produce EGSnrc transport text, for inclusion in EGSnrc input files.

    Args:
        electron_energy_cutoff (float): The global ECUT value.
        photon_energy_cutoff (float): The global PCUT value.
        exact_boundary_crossing (bool): Flag for EXACT boundary crossing algorithm use.
        low_energy_simulation (bool): Flag for simulation of low energy particles.
    """
    parameter_text = []
    parameter_text.append(" #########################")
    parameter_text.append(" :Start MC Transport Parameter:")
    parameter_text.append("")
    parameter_text.append(" Global ECUT= {:.5}".format(electron_energy_cutoff))
    parameter_text.append(" Global PCUT= {:.5}".format(photon_energy_cutoff))
    if exact_boundary_crossing:
        parameter_text.append(" Global SMAX= 1E10")
    else:
        parameter_text.append(" Global SMAX= 5")
    parameter_text.append(" ESTEPE= 0.25")
    parameter_text.append(" XIMAX= 0.5")
    if exact_boundary_crossing:
        parameter_text.append("Boundary crossing algorithm= EXACT")
        parameter_text.append("Skin depth for BCA= 3")
    else:
        parameter_text.append("Boundary crossing algorithm= PRESTA-I")
        parameter_text.append("Skin depth for BCA= 0")
    parameter_text.append(" Electron-step algorithm= PRESTA-II");
    parameter_text.append(" Spin effects= On")
    if low_energy_simulation:
        parameter_text.append(" Brems angular sampling= KM")
    else:
        parameter_text.append(" Brems angular sampling= Simple");
    parameter_text.append(" Brems cross sections= BH")
    parameter_text.append(" Triplet production= Off")
    parameter_text.append(" Bound Compton scattering= On")
    parameter_text.append(" Radiative Compton corrections= Off")
    parameter_text.append(" Pair angular sampling= Simple")
    if low_energy_simulation:
        parameter_text.append(" Pair cross sections= NRC");
    else:
        parameter_text.append(" Pair cross sections= BH")
    parameter_text.append(" Photoelectron angular sampling= Off");
    if low_energy_simulation:
        parameter_text.append(" Rayleigh scattering= On");
    else:
        parameter_text.append(" Rayleigh scattering= Off");
    if low_energy_simulation:
        parameter_text.append(" Atomic relaxations= On");
    else:
        parameter_text.append(" Atomic relaxations= Off");
    parameter_text.append(" Electron impact ionization= Off")
    parameter_text.append("")
    parameter_text.append(" :Stop MC Transport Parameter:")
    parameter_text.append(" #########################")
    return parameter_text
