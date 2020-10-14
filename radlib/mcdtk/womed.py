# -*- coding: utf-8 -*-
""" WOmed T-300 EGSnrc file writing module.

This module provides BEAMnrc input file writing functionality for simulations of the a WOmed T-300 system.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = []
__license__ = 'GPL3'

# import required code
import os

# define default materials
_default_air_material = 'AIR521ICRU'
_default_al_material = 'AL521ICRU'
_default_be_material = 'BE521ICRU'
_default_cu_material = 'CU521ICRU'
_default_h2o_material = 'H2O521ICRU'
_default_pb_material = 'PB521ICRU'
_default_sn_material = 'SN521ICRU'
_default_w_material = 'W521ICRU'

# define detault EGSnrc parameters
_default_ecut = '0.521'
_default_pcut = '0.01'
_default_bca = 'EXACT'
_default_rayleigh_flag = 'On'
_default_eii_flag = 'On'
_default_bcse_enhancement = 100

# define detault job running parameters
_default_pegs_data = '521womed'
_default_queue = 'pbs_medium'
_default_jobs = 10

# define available components
_filters = [
    { _default_al_material : 0.06 },
    { _default_al_material : 0.18 },
    { _default_al_material : 0.34 },
    { _default_cu_material : 0.02, _default_al_material : 0.05 },
    { _default_cu_material : 0.1, _default_al_material : 0.1 },
    { _default_cu_material : 0.16, _default_al_material : 0.05 },
    { _default_sn_material : 0.05, _default_cu_material : 0.26, _default_al_material : 0.05}
]
_circular_applicators = { 'd2' : [1, 30], 'd3' : [1.5, 30], 'd5' : [2.5, 30], 'd10' : [5, 30]}
_circular_endplates = { } 
_square_applicators = { '5x7' : [5, 7, 50], '8x8' : [8, 8, 50], '10x10' : [10, 10, 50], '10x20' : [10, 20, 50], '15x15' : [15, 15, 50], '20x20' : [20, 20, 50]}
_square_endplates = { }

def write_beamnrc_input(path : str, simulation_name : str, energy : float, beam_filter : int, applicator: str, histories=1000000000):
    """ Creates BEAMnrc input file for a WOmed T-300 simulation.

    Args:
       path (str): The path of the input file to be created.
       simulation_name (str): The name of the simulation, to be inserted into BEAMnrc input file.
       energy (float): The energy of the electron beam incident on target, in MeV.
       beam_filter (int): The filter to use.
       applicator (str): The applicator to use.
       histories (int): The number of histories to use in the simulation.
    """
    # pre-conditions
    if beam_filter < 1 or beam_filter > len(_filters):
        raise ValueError('Invalid filter specified.')
    if applicator not in _circular_applicators and applicator not in _square_applicators:
        raise ValueError('Invalid applicator specified.')
    # open file and input instructions
    file_handle = open(path, 'w')
    file_handle.write(simulation_name + '\n')
    file_handle.write(_default_air_material + '\n')
    file_handle.write('0, 0, 0, 0, 0, 2, 0,  IWATCH ETC.\n')
    file_handle.write(str(histories) + ', 33, 97, 99, 2, 5000, 0, 0,  NCASE ETC.\n')
    file_handle.write('5, 31.5, 0, 0, 0, ,  DIRECTIONAL BREM OPTIONS\n')
    file_handle.write('-1, 10, 0.35, 1, 0, 0,  0.0, 0.0, 0.0, 0.0,  IQIN, ISOURCE + OPTIONS\n')
    file_handle.write('0, MONOENERGETIC\n')
    file_handle.write(str(energy) + '\n')
    file_handle.write('0, 0, ' + _default_ecut + ', ' + _default_pcut + ', 0, 0, ,  0 , ECUT,PCUT,IREJCT,ESAVE\n')
    file_handle.write('0, 0, 0, 0, 0,  PHOTON FORCING\n')
    file_handle.write('1, 8,  SCORING INPUT\n')
    file_handle.write('1, 1\n')
    file_handle.write('3, \n')
    file_handle.write('0,  DOSE COMPONENTS\n')
    file_handle.write('0.0, Z TO FRONT FACE\n')
    # write beam production information
    file_handle.write('*********** start of CM XTUBE with identifier TUBE  ***********\n')
    file_handle.write('15, RMAX\n')
    file_handle.write('TUBE\n')
    file_handle.write('0, 3, ZMIN, ZTHICK\n')
    file_handle.write('30, ANGLE\n')
    file_handle.write('1, # LAYERS\n')
    file_handle.write('0.3, 0\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0, \n')
    file_handle.write(_default_w_material + '\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0, \n')
    file_handle.write('VACUUM\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0, \n')
    file_handle.write(_default_cu_material + '\n')
    file_handle.write('*********** start of CM SLABS with identifier VACUUM  ***********\n')
    file_handle.write('15, RMAX\n')
    file_handle.write('VACUUM\n')
    file_handle.write('1, NSLABS\n')
    file_handle.write('3, ZMIN\n')
    file_handle.write('3.5, ' + _default_ecut + ', ' + _default_pcut + ', 0, 0, 0\n')
    file_handle.write('VACUUM\n')
    file_handle.write('*********** start of CM SLABS with identifier XWIN  ***********\n')
    file_handle.write('15, RMAX\n')
    file_handle.write('XWIN\n')
    file_handle.write('1, NSLABS\n')
    file_handle.write('6.5, ZMIN\n')
    file_handle.write('0.3, ' + _default_ecut + ', ' + _default_pcut + ', 0, 0, 0\n')
    file_handle.write(_default_be_material + '\n')
    file_handle.write('*********** start of CM CONS3R with identifier PRIMCOLL  ***********\n')
    file_handle.write('15, RMAX\n')
    file_handle.write('PRIMCOLL\n')
    file_handle.write('6.8, ZMIN\n')
    file_handle.write('2.3, ZTHICK\n')
    file_handle.write('6, NUM_NODE\n')
    file_handle.write('6.8, 2, \n')
    file_handle.write('7.0, 2, \n')
    file_handle.write('7.8, 2.3, \n')
    file_handle.write('8.8, 3.25, \n')
    file_handle.write('8.8, 3.65, \n')
    file_handle.write('9.1, 3.65, \n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0, 0, \n')
    file_handle.write(_default_air_material + '\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0, 0, \n')
    file_handle.write(_default_pb_material + '\n')
    # write filter information
    file_handle.write('*********** start of CM SLABS with identifier FILTER  ***********\n')
    file_handle.write('15, RMAX\n')
    file_handle.write('FILTER\n')
    num_filter_materials = len(_filters[beam_filter - 1])
    file_handle.write(str(num_filter_materials) + ', NSLABS\n')
    file_handle.write('12, ZMIN\n')
    for filter_material in _filters[beam_filter - 1]:
        filter_material_thickness = _filters[beam_filter - 1][filter_material]
        file_handle.write(str(filter_material_thickness) + ', ' + _default_ecut + ', ' + _default_pcut + ', 0, , 0\n')
        file_handle.write(filter_material + '\n')
    # write monitor chamber information
    file_handle.write('*********** start of CM SLABS with identifier MONCHAM  ***********\n')
    file_handle.write('15, RMAX\n')
    file_handle.write('MONCHAM\n')
    file_handle.write('1, NSLABS\n')
    file_handle.write('13, ZMIN\n')
    file_handle.write('0.03, ' + _default_ecut + ', ' + _default_pcut + ', 0, , 0\n')
    file_handle.write(_default_al_material + '\n')
    # write applicator information
    if applicator in _circular_applicators:
        field_radius = _circular_applicators[applicator][0]
        focal_spot_distance = _circular_applicators[applicator][1]
        file_handle.write('*********** start of CM CONS3R with identifier APERT  ***********\n')
        file_handle.write('15, RMAX\n')
        file_handle.write('APERT\n')
        file_handle.write('16.7, ZMIN\n')
        file_handle.write('0.8, ZTHICK\n')
        file_handle.write('2, NUM_NODE\n')
        file_handle.write('16.7, ' + str(round(field_radius * (16.7-1.5) / focal_spot_distance, 5)) + ', \n')
        file_handle.write('17.5, ' + str(round(field_radius * (17.5-1.5) / focal_spot_distance, 5)) + ', \n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n')
        file_handle.write(_default_air_material + '\n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n')
        file_handle.write(_default_pb_material + '\n')
        file_handle.write('*********** start of CM CONS3R with identifier APPL ***********\n')
        file_handle.write('3, RMAX\n')
        file_handle.write('APPL\n')
        file_handle.write('17.5, ZMIN\n')
        file_handle.write(str(focal_spot_distance - 17.5 + 1.5) + ', ZTHICK\n')
        file_handle.write('2, NUM_NODE\n')
        file_handle.write('17.5, ' + str(round(field_radius, 5)) + ', \n')
        file_handle.write(str(focal_spot_distance + 1.5) + ', ' + str(round(field_radius, 5)) + ', \n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n')
        file_handle.write(_default_air_material + '\n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n')
        file_handle.write(_default_pb_material + '\n')
        # add end plate as necessary
        if applicator in _circular_endplates:
            # TODO add end plate
            print('Not yet implemented')
    # if rectangular applicator
    elif applicator in _square_applicators:
        # TODO model the square applicators
        print('Not yet implemented')
        if applicator in _square_endplates:
            # TODO add end plate
            print('Not yet implemented')
    file_handle.write('*********************end of all CMs*****************************\n')
    file_handle.write('#########################\n')
    file_handle.write(':Start MC Transport Parameter:\n')
    file_handle.write('    \n')
    file_handle.write('Global ECUT= ' + _default_ecut + '\n')
    file_handle.write('Global PCUT= ' + _default_pcut + '\n')
    file_handle.write('Global SMAX= 1e10\n')
    file_handle.write('ESTEPE= 0.25\n')
    file_handle.write('XIMAX= 0.5\n')
    file_handle.write('Boundary crossing algorithm= ' + _default_bca + '\n')
    file_handle.write('Skin depth for BCA= 0\n')
    file_handle.write('Electron-step algorithm= PRESTA-II\n')
    file_handle.write('Spin effects= On\n')
    file_handle.write('Brems angular sampling= Simple\n')
    file_handle.write('Brems cross sections= BH\n')
    file_handle.write('Bound Compton scattering= Norej\n')
    file_handle.write('Compton cross sections= default\n')
    file_handle.write('Pair angular sampling= Simple\n')
    file_handle.write('Pair cross sections= BH\n')
    file_handle.write('Photoelectron angular sampling= Off\n')
    file_handle.write('Rayleigh scattering= ' + _default_rayleigh_flag + '\n')
    file_handle.write('Atomic relaxations= On\n')
    file_handle.write('Electron impact ionization= ' + _default_eii_flag + '\n')
    file_handle.write('Photon cross sections= PEGS4\n')
    file_handle.write('Photon cross-sections output= Off\n')
    file_handle.write('    \n')
    file_handle.write(':Stop MC Transport Parameter:\n')
    file_handle.write('#########################\n')
    file_handle.write(':Start BCSE:\n')
    file_handle.write('    \n')
    file_handle.write('Use BCSE= On\n')
    file_handle.write('Media to enhance=  ' + _default_w_material + '\n')
    file_handle.write('Enhancement constant= ' + str(_default_bcse_enhancement) + '\n')
    file_handle.write('Enhancement power= 0\n')
    file_handle.write('Russian Roulette= off\n')
    file_handle.write('\n')
    file_handle.write(':Stop BCSE:\n')
    file_handle.write('#########################\n')