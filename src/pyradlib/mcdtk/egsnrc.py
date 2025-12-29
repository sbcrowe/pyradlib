# -*- coding: utf-8 -*-
"""EGSnrc module.

This module provides EGSnrc file handling functionality.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import numpy as np
from numpy.typing import ArrayLike


def read_egsphant(path: str):
    """Read content of EGSnrc EGSPHANT files.

    Parameters
    ----------
    path : str
        Path of file containing EGSPHANT data.

    Returns
    -------
    material_names : dict
        Dictionary containing material names and associated media assignment numbers.
    x_boundaries, y_boundaries, z_boundaries : list
        Phantom geometry voxel boundary coordinates.
    media, densities : array_like
        Phantom geometry material assignments and density values.
    """
    with open(path) as f:
        num_materials = int(f.readline().strip())
        material_names = {}
        for i in range(num_materials):
            material_names[i + 1] = f.readline().strip()
        f.readline()
        num_voxels = [int(x) for x in f.readline().strip().split()]
        x_boundaries = [float(x) for x in f.readline().strip().split()]
        y_boundaries = [float(y) for y in f.readline().strip().split()]
        z_boundaries = [float(z) for z in f.readline().strip().split()]
        media = np.zeros((num_voxels[2], num_voxels[1], num_voxels[0]))
        for k in range(num_voxels[2]):
            for j in range(num_voxels[1]):
                media[k, j] = [int(x) for x in f.readline().strip()]
            f.readline()
        densities = np.zeros((num_voxels[2], num_voxels[1], num_voxels[0]))
        for k in range(num_voxels[2]):
            for j in range(num_voxels[1]):
                densities[k, j] = [float(x) for x in f.readline().strip().split()]
            f.readline()
        return (
            material_names,
            x_boundaries,
            y_boundaries,
            z_boundaries,
            media,
            densities,
        )


# egsphant file write
def write_egsphant(
    path: str,
    material_names: dict[int, str],
    x_boundaries: list[float],
    y_boundaries: list[float],
    z_boundaries: list[float],
    media: ArrayLike,
    densities: ArrayLike,
):
    with open(path, "w") as f:
        num_materials = len(material_names)
        f.write(str(num_materials) + "\n")
        for material_name in list(material_names.values()):
            f.write(material_name + "\n")
        f.write(" ".join(["1.00000000"] * 6) + "\n")
        num_x = len(x_boundaries) - 1
        num_y = len(y_boundaries) - 1
        num_z = len(z_boundaries) - 1
        f.write(" " + str(num_x) + " " + str(num_y) + " " + str(num_z) + "\n")
        f.write(" ".join([str(x) for x in x_boundaries]) + "\n")
        f.write(" ".join([str(y) for y in y_boundaries]) + "\n")
        f.write(" ".join([str(z) for z in z_boundaries]) + "\n")
        for k in range(num_z):
            for j in range(num_y):
                f.write("".join([str(int(x)) for x in media[k, j]]) + "\n")
            f.write("\n")
        for k in range(num_z):
            for j in range(num_y):
                f.write(" ".join([str(x) for x in densities[k, j]]) + "\n")
            f.write("\n")


def transport_parameter_text(
    electron_energy_cutoff: float = 0.521,
    photon_energy_cutoff: float = 0.01,
    exact_boundary_crossing: bool = True,
    low_energy_simulation: bool = False,
) -> list[str]:
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
