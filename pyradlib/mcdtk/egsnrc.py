# -*- coding: utf-8 -*-
""" EGSnrc module.

This module provides EGSnrc file handling functionality.
"""


def transport_parameter_text(
    electron_energy_cutoff=0.521,
    photon_energy_cutoff=0.01,
    exact_boundary_crossing=True,
    low_energy_simulation=False,
):
    """Produce EGSnrc transport text, for inclusion in EGSnrc input files.

    Parameters
    ----------
    electron_energy_cutoff : float, optional
        Global ECUT value.
    photon_energy_cutoff : float
        Global PCUT value.
    exact_boundary_crossing : bool
        Flag for use of EXACT boundary crossing algorithm use (default is True).
    low_energy_simulation : bool
        Flag for simulation of low energy particles, including NRC cross sections, Rayleigh scattering, atomic relaxation, etc. (default is False).

    Returns
    -------
    transport_parameters : array_like
        List containing EGSnrc transport parameter strings, for us in EGSnrc input file.
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
    parameter_text.append(" Electron-step algorithm= PRESTA-II")
    parameter_text.append(" Spin effects= On")
    if low_energy_simulation:
        parameter_text.append(" Brems angular sampling= KM")
    else:
        parameter_text.append(" Brems angular sampling= Simple")
    parameter_text.append(" Brems cross sections= BH")
    parameter_text.append(" Triplet production= Off")
    parameter_text.append(" Bound Compton scattering= On")
    parameter_text.append(" Radiative Compton corrections= Off")
    parameter_text.append(" Pair angular sampling= Simple")
    if low_energy_simulation:
        parameter_text.append(" Pair cross sections= NRC")
    else:
        parameter_text.append(" Pair cross sections= BH")
    parameter_text.append(" Photoelectron angular sampling= Off")
    if low_energy_simulation:
        parameter_text.append(" Rayleigh scattering= On")
    else:
        parameter_text.append(" Rayleigh scattering= Off")
    if low_energy_simulation:
        parameter_text.append(" Atomic relaxations= On")
    else:
        parameter_text.append(" Atomic relaxations= Off")
    parameter_text.append(" Electron impact ionization= Off")
    parameter_text.append("")
    parameter_text.append(" :Stop MC Transport Parameter:")
    parameter_text.append(" #########################")
    return parameter_text
