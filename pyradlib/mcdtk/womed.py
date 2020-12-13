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
import random
import numpy as np

# define default materials
_default_air_material = 'AIR521ICRU'
_default_al_material = 'AL521ICRU'
_default_be_material = 'BE521ICRU'
_default_cu_material = 'CU521ICRU'
_default_h2o_material = 'H2O521ICRU'
_default_pb_material = 'PB521ICRU'
_default_pmma_material = 'PMMA521ICRU'
_default_sn_material = 'SN521ICRU'
_default_w_material = 'W521ICRU'

# define detault EGSnrc parameters
_default_ecut = '0.521'
_default_pcut = '0.01'
_default_bca = 'EXACT'
_default_rayleigh_flag = 'On'
_default_eii_flag = 'On'
_default_bcse_enhancement = 100
_default_beamnrc_histories = 1000000000
_default_dosxyznrc_histories = 25000000000

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
_cylindrical_applicators = { 'd2' : [1, 30], 'd3' : [1.5, 30], 'd5' : [2.5, 30] }
_cylindrical_endplates = { } 
_conical_applicators = { 'd10' : [5, 30] }
_conical_endplates = { 'd10' : 0.1 }
_square_applicators = { '5x7' : [5, 7, 50], '8x8' : [8, 8, 50], '10x10' : [10, 10, 50], '10x20' : [10, 20, 50], '15x15' : [15, 15, 50], '20x20' : [20, 20, 50]}
_square_endplates = { '5x7' : 0.1, '8x8' : 0.1, '10x10' : 0.1, '10x20' : 0.1, '15x15' : 0.1, '20x20' : 0.1 }
_applicator_wall = 0.5

def write_beamnrc_input(path : str, simulation_name : str, energy : float, beam_filter : int, applicator: str, histories=_default_beamnrc_histories):
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
    if applicator not in _cylindrical_applicators and applicator not in _conical_applicators and applicator not in _square_applicators:
        raise ValueError('Invalid applicator specified.')
    # open file and input instructions
    file_handle = open(path, 'w')
    file_handle.write(simulation_name + '\n')
    file_handle.write(_default_air_material + '\n')
    file_handle.write('0, 0, 0, 0, 0, 2, 0,  IWATCH ETC.\n')
    ixxin = random.randint(0,4)
    jxxin = random.randint(0,30081)
    file_handle.write(str(histories) + ', ' + str(ixxin) + ', ' + str(jxxin) + ', 99, 2, 5000, 0, 0,  NCASE ETC.\n')
    if applicator in _cylindrical_applicators:
        dbs_fs = _cylindrical_applicators[applicator][0]
        dbs_ssd = _cylindrical_applicators[applicator][1] + 1.5
    elif applicator in _conical_applicators:
        dbs_fs = _conical_applicators[applicator][0]
        dbs_ssd = _conical_applicators[applicator][1] + 1.5
    elif applicator in _square_applicators:
        dbs_fs = np.max(_square_applicators[applicator][0:-1])
        dbs_ssd = _square_applicators[applicator][2]
    file_handle.write(str(dbs_fs) + ', ' + str(dbs_ssd) + ', 0, 0, 0, ,  DIRECTIONAL BREM OPTIONS\n')
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
    if applicator in _cylindrical_applicators:
        field_radius = _cylindrical_applicators[applicator][0]
        focal_spot_distance = _cylindrical_applicators[applicator][1]
        if applicator in _cylindrical_endplates:
            end_plate_thickness = _cylindrical_endplates[applicator]
        else:
            end_plate_thickness = 0
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
        file_handle.write(str(round(field_radius + _applicator_wall, 5)) + ', RMAX\n')
        file_handle.write('APPL\n')
        file_handle.write('17.5, ZMIN\n')
        file_handle.write(str(focal_spot_distance - 17.5 + 1.5 - end_plate_thickness) + ', ZTHICK\n')
        file_handle.write('2, NUM_NODE\n')
        file_handle.write('17.5, ' + str(round(field_radius, 5)) + ', \n')
        file_handle.write(str(focal_spot_distance + 1.5 - end_plate_thickness) + ', ' + str(round(field_radius, 5)) + ', \n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n')
        file_handle.write(_default_air_material + '\n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n')
        file_handle.write(_default_pb_material + '\n')
        # add end plate as necessary
        if applicator in _cylindrical_endplates:
            file_handle.write('*********** start of CM SLABS with identifier PLATE ***********\n')
            file_handle.write(str(round(field_radius + _applicator_wall, 5)) + ', RMAX\n')
            file_handle.write('PLATE\n')
            file_handle.write('1, NSLABS\n')
            file_handle.write(str(focal_spot_distance + 1.5 - end_plate_thickness) + ', ZMIN\n')
            file_handle.write(str(end_plate_thickness) + ', ' + _default_ecut + ', ' + _default_pcut + ', 0, , 0\n')
            file_handle.write(_default_pmma_material + '\n')
    # if conical applicator
    elif applicator in _conical_applicators:
        field_radius = _conical_applicators[applicator][0]
        focal_spot_distance = _conical_applicators[applicator][1]
        if applicator in _conical_endplates:
            end_plate_thickness = _conical_endplates[applicator]
        else:
            end_plate_thickness = 0
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
        file_handle.write('*********** start of CM FLATFILT with identifier APPL  ***********\n')
        file_handle.write('15, RMAX\n') # RMAX_CM(ICM_$FLATFILT) (F10.0): Radius of outer boundary of CM (cm).
        file_handle.write('APPL\n')
        file_handle.write('17.5, ZMIN\n')
        file_handle.write('1,\n') # no. layers
        file_handle.write('2,' + str(focal_spot_distance - 17.5 + 1.5 - end_plate_thickness) + ',\n') # no. cones, thickness
        file_handle.write(str(round(field_radius * 0.9, 5)) + ',' + str(round(field_radius * 0.9 + _applicator_wall, 5)) + ',\n') # rtop
        file_handle.write(str(round(field_radius, 5)) + ',' + str(round(field_radius + _applicator_wall, 5)) + ',\n') # rbot
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # cone 1
        file_handle.write(_default_air_material + '\n') # cone 1
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # cone 2
        file_handle.write(_default_pb_material + '\n') # cone 2
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # outside cone
        file_handle.write(_default_air_material + '\n') # outside cone
        if applicator in _conical_endplates:
            file_handle.write('*********** start of CM SLABS with identifier PLATE ***********\n')
            file_handle.write(str(round(field_radius + _applicator_wall, 5)) + ', RMAX\n')
            file_handle.write('PLATE\n')
            file_handle.write('1, NSLABS\n')
            file_handle.write(str(focal_spot_distance + 1.5 - end_plate_thickness) + ', ZMIN\n')
            file_handle.write(str(end_plate_thickness) + ', ' + _default_ecut + ', ' + _default_pcut + ', 0, , 0\n')
            file_handle.write(_default_pmma_material + '\n')
    # if rectangular applicator
    elif applicator in _square_applicators:
        field_length_1 = _square_applicators[applicator][0]
        field_length_2 = _square_applicators[applicator][1]
        focal_spot_distance = _square_applicators[applicator][2]
        if applicator in _conical_endplates:
            end_plate_thickness = _conical_endplates[applicator]
        else:
            end_plate_thickness = 0
        file_handle.write('*********** start of CM PYRAMIDS with identifier APERT  ***********\n')
        file_handle.write('15, RMAX\n')
        file_handle.write('APERT\n')
        file_handle.write('1,0\n')
        apert_xf = round((field_length_1 / 2) * (16.7-1.5) / focal_spot_distance, 5)
        apert_xb = round((field_length_1 / 2) * (17.5-1.5) / focal_spot_distance, 5)
        apert_yf = round((field_length_2 / 2) * (16.7-1.5) / focal_spot_distance, 5)
        apert_yb = round((field_length_2 / 2) * (17.5-1.5) / focal_spot_distance, 5)
        file_handle.write('16.7, 17.5, ' + str(apert_xf) + ', ' + str(apert_xb) + ', ' + str(-apert_xf) + ', ' + str(-apert_xb) + ', ' + str(apert_yf) + ', ' + str(apert_yb) + ', ' + str(-apert_yf) + ', ' + str(-apert_yb) + ', 15, 15\n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # air
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # flange
        file_handle.write(_default_pb_material + '\n')
        file_handle.write('*********** start of CM PYRAMIDS with identifier APPL  ***********\n')
        file_handle.write('15, RMAX\n')
        file_handle.write('APPL\n')
        file_handle.write('1,0\n')
        appl_xf = round((field_length_1 / 2) * (17.6-1.5) / focal_spot_distance, 5)
        appl_xb = round((field_length_1 / 2) * (focal_spot_distance - end_plate_thickness) / focal_spot_distance, 5)
        appl_yf = round((field_length_2 / 2) * (17.6-1.5) / focal_spot_distance, 5)
        appl_yb = round((field_length_2 / 2) * (focal_spot_distance - end_plate_thickness) / focal_spot_distance, 5)
        file_handle.write('17.6, ' + str(focal_spot_distance - end_plate_thickness) + ', ' + str(appl_xf) + ', ' + str(appl_xb) + ', ' + str(-appl_xf) + ', ' + str(-appl_xb) + ', ' + str(appl_yf) + ', ' + str(appl_yb) + ', ' + str(-appl_yf) + ', ' + str(-appl_yb) + ', ' + str((field_length_1 / 2) + _applicator_wall) + ', ' + str((field_length_2 / 2) + _applicator_wall) + '\n')
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # air
        file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, , 0, \n') # applicator
        file_handle.write(_default_pb_material + '\n')
        if applicator in _square_endplates:
            file_handle.write('*********** start of CM SLABS with identifier PLATE ***********\n')
            file_handle.write(str(round(max(field_length_1 / 2, field_length_2 / 2) + _applicator_wall, 5)) + ', RMAX\n')
            file_handle.write('PLATE\n')
            file_handle.write('1, NSLABS\n')
            file_handle.write(str(focal_spot_distance + 1.5 - end_plate_thickness) + ', ZMIN\n')
            file_handle.write(str(end_plate_thickness) + ', ' + _default_ecut + ', ' + _default_pcut + ', 0, , 0\n')
            file_handle.write(_default_pmma_material + '\n')
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

def write_dosxyznrc_pdd_input(path : str, phase_space_path : str, applicator : str, simulation_name : str = 'PDD simulation', histories=_default_dosxyznrc_histories):
    """ Creates DOSXYZnrc input file for PDD simulations in a 40x40x40 water phantom.

    Args:
        path (str): The path of the input file to be created.
        phase_space_path (str): The path of the phase space file to be used.
        simulation_name (str): The name of the simulation, to be inserted into DOSXYZnrc input file.
        histories (int): The number of histories to use in the simulation.
    """
    file_handle = open(path, 'w')
    file_handle.write(simulation_name + '\n')
    file_handle.write('1\n')
    file_handle.write(_default_h2o_material + '\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0\n')
    file_handle.write('3, 3, -4, 0\n')
    file_handle.write('-20\n')
    file_handle.write('-0.1\n')
    file_handle.write('0.1\n')
    file_handle.write('20\n')
    file_handle.write('-20\n')
    file_handle.write('-0.1\n')
    file_handle.write('0.1\n')
    file_handle.write('20\n')
    file_handle.write('0.0\n')
    file_handle.write('0.01, 15\n')
    file_handle.write('0.1, 99\n')
    file_handle.write('1, 10\n')
    file_handle.write('10, 2\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    file_handle.write('2, 2, 2, 2, 0, 200, 1, 0\n') # egslst printing
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    if applicator in _cylindrical_applicators:
        dbs_fs = _cylindrical_applicators[applicator][0]
        dbs_ssd = _cylindrical_applicators[applicator][1] + 1.5
    elif applicator in _conical_applicators:
        dbs_fs = _conical_applicators[applicator][0]
        dbs_ssd = _conical_applicators[applicator][1] + 1.5
    elif applicator in _square_applicators:
        dbs_fs = np.max(_square_applicators[applicator][0:-1])
        dbs_ssd = _square_applicators[applicator][2]
    file_handle.write('2, 2, 0, 0, 0, 180, 90, 0, 180, 1, ' + str(dbs_fs) + ', ' + str(dbs_ssd) + ', ' + str(dbs_ssd) + ', 0\n') # iqin, etc.
    file_handle.write('2, 0, 0, 50, 0, 0, 0, 0\n')
    file_handle.write(phase_space_path + '\n')
    file_handle.write(str(histories) + ', 0, 100, 33, 97, 40, 0, 0, 0, 0, , 0, 0, 0, 1, 0, 0\n')
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
    file_handle.write('Bound Compton scattering= Off\n')
    file_handle.write('Compton cross sections= default\n')
    file_handle.write('Pair angular sampling= Simple\n')
    file_handle.write('Pair cross sections= BH\n')
    file_handle.write('Photoelectron angular sampling= Off\n')
    file_handle.write('Rayleigh scattering= ' + _default_rayleigh_flag + '\n')
    file_handle.write('Atomic relaxations= On\n')
    file_handle.write('Electron impact ionization= ' + _default_eii_flag + '\n')
    file_handle.write('Photon cross sections= PEGS4\n')
    file_handle.write('Photon cross-sections output= Off\n')
    file_handle.write('\n')
    file_handle.write(':Stop MC Transport Parameter:\n')
    file_handle.write('#########################\n')
    file_handle.write('\n')
    file_handle.close()

def write_dosxyznrc_isodose_input(path : str, phase_space_path : str, applicator : str, simulation_name : str = 'Isodose simulation', histories=_default_dosxyznrc_histories, plane = 'in-plane', phantom = 40, resolution = '0.2'):
    """ Creates DOSXYZnrc input file for in-plane isodose simulations in a 40x40x40 water phantom.

    Args:
        path (str): The path of the input file to be created.
        phase_space_path (str): The path of the phase space file to be used.
        simulation_name (str): The name of the simulation, to be inserted into DOSXYZnrc input file.
        histories (int): The number of histories to use in the simulation.
        plane (str): The dose plane to be simulated, 'in-plane' or 'cross-plane'.
        phantom (float): The size of the phantom, in cm.
        resolution (float): The resolution of the phantom, in cm.
    """
    # pre-conditions
    if plane not in ['in-plane', 'cross-plane']:
        raise ValueError('Invalid plane selection.')
    # write input file
    file_handle = open(path, 'w')
    file_handle.write(simulation_name + '\n')
    file_handle.write('1\n')
    file_handle.write(_default_h2o_material + '\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0\n')
    plane_voxels = int(phantom / resolution)
    phantom_edge = phantom / 2
    cax_edge = resolution / 2
    if plane == 'in-plane':
        file_handle.write('-1, 3, -1, 0\n')
        file_handle.write('-' + str(phantom_edge) + '\n')
        file_handle.write(str(resolution) + ', ' + str(plane_voxels) + '\n')
        file_handle.write('-' + str(phantom_edge) + '\n')
        file_handle.write('-' + str(cax_edge) + '\n')
        file_handle.write(str(cax_edge) + '\n')
        file_handle.write(str(phantom_edge) + '\n')
    else:
        file_handle.write('3, -1, -1, 0\n')
        file_handle.write('-' + str(phantom_edge) + '\n')
        file_handle.write('-' + str(cax_edge) + '\n')
        file_handle.write(str(cax_edge) + '\n')
        file_handle.write(str(phantom_edge) + '\n')
        file_handle.write('-' + str(phantom_edge) + '\n')
        file_handle.write(str(resolution) + ', ' + str(plane_voxels) + '\n')
    file_handle.write('0.0\n')
    file_handle.write(str(resolution) + ', ' + str(plane_voxels) + '\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    file_handle.write('2, 2, 2, 2, 0, 200, 1, 0\n') 
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n') 
    if applicator in _cylindrical_applicators:
        dbs_fs = _cylindrical_applicators[applicator][0]
        dbs_ssd = _cylindrical_applicators[applicator][1] + 1.5
    elif applicator in _conical_applicators:
        dbs_fs = _conical_applicators[applicator][0]
        dbs_ssd = _conical_applicators[applicator][1] + 1.5
    elif applicator in _square_applicators:
        dbs_fs = np.max(_square_applicators[applicator][0:-1])
        dbs_ssd = _square_applicators[applicator][2]
    file_handle.write('2, 2, 0, 0, 0, 180, 90, 0, 180, 1, ' + str(dbs_fs) + ', ' + str(dbs_ssd) + ', ' + str(dbs_ssd) + ', 0\n')
    file_handle.write('2, 0, 0, 50, 0, 0, 0, 0\n')
    file_handle.write(phase_space_path + '\n')
    file_handle.write(str(histories) + ', 0, 100, 33, 97, 40, 0, 0, 0, 0, , 0, 0, 0, 1, 0, 0\n')
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
    file_handle.write('Bound Compton scattering= Off\n')
    file_handle.write('Compton cross sections= default\n')
    file_handle.write('Pair angular sampling= Simple\n')
    file_handle.write('Pair cross sections= BH\n')
    file_handle.write('Photoelectron angular sampling= Off\n')
    file_handle.write('Rayleigh scattering= ' + _default_rayleigh_flag + '\n')
    file_handle.write('Atomic relaxations= On\n')
    file_handle.write('Electron impact ionization= ' + _default_eii_flag + '\n')
    file_handle.write('Photon cross sections= PEGS4\n')
    file_handle.write('Photon cross-sections output= Off\n')
    file_handle.write('\n')
    file_handle.write(':Stop MC Transport Parameter:\n')
    file_handle.write('#########################\n')
    file_handle.write('\n')
    file_handle.close()

def write_uniform_dose_input(path : str, voxel_size : float, phantom_size : float, phase_space_path : str, applicator : str, simulation_name : str = 'Dose simulation', histories=_default_dosxyznrc_histories):
    """ Creates DOSXYZnrc input file for dose simulations in a uniform resolution cubic water phantom.

    Args:
        path (str): The path of the input file to be created.
        voxel_size (float): The size of the voxels, in cm.
        phantom_size (float): The size of the phantom, in cm.
        phase_space_path (str): The path of the phase space file to be used.
        simulation_name (str): The name of the simulation, to be inserted into DOSXYZnrc input file.
        histories (int): The number of histories to use in the simulation.
    """
    # pre-conditions
    if not np.isclose(phantom_size % voxel_size, 0):
        raise ValueError('Phantom size not divisible by voxel size.')
    # determine number of voxels
    boundary = phantom_size / 2
    voxels = int(phantom_size / voxel_size)
    # write file
    file_handle = open(path, 'w')
    file_handle.write(simulation_name + '\n')
    file_handle.write('1\n')
    file_handle.write(_default_h2o_material + '\n')
    file_handle.write(_default_ecut + ', ' + _default_pcut + ', 0, 0\n')
    file_handle.write('-1, -1, -1, 0\n')
    file_handle.write('-' + boundary + '\n')
    file_handle.write(str(voxel_size) + ', ' + str(voxels) + '\n')
    file_handle.write('-' + boundary + '\n')
    file_handle.write(str(voxel_size) + ', ' + str(voxels) + '\n')
    file_handle.write('0.0\n')
    file_handle.write(str(voxel_size) + ', ' + str(voxels) + '\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    file_handle.write('2, 2, 2, 2, 0, 200, 1, 0\n')
    file_handle.write('0, 0, 0, 0, 0, 0, 0, 0\n')
    if applicator in _cylindrical_applicators:
        dbs_fs = _cylindrical_applicators[applicator][0]
        dbs_ssd = _cylindrical_applicators[applicator][1] + 1.5
    elif applicator in _conical_applicators:
        dbs_fs = _conical_applicators[applicator][0]
        dbs_ssd = _conical_applicators[applicator][1] + 1.5
    elif applicator in _square_applicators:
        dbs_fs = np.max(_square_applicators[applicator][0:-1])
        dbs_ssd = _square_applicators[applicator][2]
    file_handle.write('2, 2, 0, 0, 0, 180, 90, 0, 180, 1, ' + str(dbs_fs) + ', ' + str(dbs_ssd) + ', ' + str(dbs_ssd) + ', 0\n')
    file_handle.write('2, 0, 0, 50, 0, 0, 0, 0\n')
    file_handle.write(phase_space_path + '\n')
    file_handle.write(str(histories) + ', 0, 100, 33, 97, 40, 0, 0, 0, 0, , 0, 0, 0, 1, 0, 0\n')
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
    file_handle.write('Bound Compton scattering= Off\n')
    file_handle.write('Compton cross sections= default\n')
    file_handle.write('Pair angular sampling= Simple\n')
    file_handle.write('Pair cross sections= BH\n')
    file_handle.write('Photoelectron angular sampling= Off\n')
    file_handle.write('Rayleigh scattering= ' + _default_rayleigh_flag + '\n')
    file_handle.write('Atomic relaxations= On\n')
    file_handle.write('Electron impact ionization= ' + _default_eii_flag + '\n')
    file_handle.write('Photon cross sections= PEGS4\n')
    file_handle.write('Photon cross-sections output= Off\n')
    file_handle.write('\n')
    file_handle.write(':Stop MC Transport Parameter:\n')
    file_handle.write('#########################\n')
    file_handle.write('\n')
    file_handle.close()