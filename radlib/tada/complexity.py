# -*- coding: utf-8 -*-
""" Treatment complexity assessment module..

This module provides radiotherapy treatment complexity assessment functionality.
"""

# authorship information
__author__ = 'Scott Crowe'
__email__ = 'sb.crowe@gmail.com'
__credits__ = ['Tanya Kairn']
__license__ = "GPL3"

# import required code
import dicom
import numpy as np
from PIL import Image
from scipy.integrate import simps

# define default parameters
_default_integration_limit = 1
_default_integration_epsilon = 0.05
_default_identifies = ['name', 'mu']
_default_metrics = ['MCS']


def aperture_area_variability(beam):
    """ Calculate the aperture area variability of a treatment beam.

    Args:
        beam (Dataset): The beam for the aperture area variability to be calculated.
    """
    result = 0
    previous_cumulative_weight = 0.0
    total_cumulative_weight = beam.FinalCumulativeMetersetWeight
    leaf_positions = []
    orthogonal_jaw_positions = []
    leaf_boundaries = _leaf_position_boundaries(beam)
    maximum_leaf_separation = _maximum_leaf_separations(beam)
    for control_point in beam.ControlPoints:
        differential_weight = (control_point.CumulativeMetersetWeight - previous_cumulative_weight) / \
                              total_cumulative_weight
        previous_cumulative_weight = control_point.CumulativeMetersetWeight
        control_point_leaf_positions = _control_point_leaf_positions(control_point)
        control_point_orthogonal_jaw_positions = _control_point_orthogonal_jaw_positions(control_point)
        if control_point_leaf_positions is not None:
            leaf_positions = control_point_leaf_positions
        if control_point_orthogonal_jaw_positions is not None:
            orthogonal_jaw_positions = control_point_orthogonal_jaw_positions
        if differential_weight > 0:
            result += _control_point_aperture_area_variability(leaf_positions, orthogonal_jaw_positions,
                                                               leaf_boundaries, maximum_leaf_separation) * \
                      differential_weight
    return result


def closed_leaf_score(beam):
    """ Calculate the closed leaf scores of a treatment beam.

    Args:
        beam (Dataset): The beam for the closed leaf score to be calculated.
    """
    return _weighted_metric(beam, 'CLS')


def cross_axis_score(beam):
    """ Calculate the cross axis scores of a treatment beam.

    Args:
        beam (Dataset): The beam for the cross axis score to be calculated.
    """
    return _weighted_metric(beam, 'CAS')


def leaf_sequence_variability(beam):
    """ Calculate the leaf sequence variability of a treatment beam.

        Args:
            beam (Dataset): The beam for the leaf sequence variability to be calculated.
        """
    return _weighted_metric(beam, 'LSV')


def mean_asymmetry_distance(beam):
    """ Calculate the mean asymmetry distance of a treatment beam.

    Args:
        beam (Dataset): The beam for the mean asymmetry distance to be calculated.
    """
    return _weighted_metric(beam, 'MAD')


def modulation_complexity_score(beam):
    """ Calculate the modulation complexity score of a treatment beam.

    Args:
        beam (Dataset): The beam for the modulation complexity score to be calculated.
    """
    result = 0
    previous_cumulative_weight = 0.0
    total_cumulative_weight = beam.FinalCumulativeMetersetWeight
    maximum_aperture_separation = _maximum_leaf_separations(beam)
    leaf_positions = []
    orthogonal_jaw_positions = []
    leaf_position_boundaries = _leaf_position_boundaries(beam)
    for control_point in beam.ControlPoints:
        differential_weight = (control_point.CumulativeMetersetWeight - previous_cumulative_weight) / \
                              total_cumulative_weight
        previous_cumulative_weight = control_point.CumulativeMetersetWeight
        control_point_leaf_positions = _control_point_leaf_positions(control_point)
        control_point_orthogonal_jaw_positions = _control_point_orthogonal_jaw_positions(control_point)
        if control_point_leaf_positions is not None:
            leaf_positions = control_point_leaf_positions
        if control_point_orthogonal_jaw_positions is not None:
            orthogonal_jaw_positions = control_point_orthogonal_jaw_positions
        if differential_weight > 0:
            control_point_aav = _control_point_aperture_area_variability(leaf_positions, orthogonal_jaw_positions,
                                                                         leaf_position_boundaries,
                                                                         maximum_aperture_separation)
            control_point_lsv = _control_point_leaf_sequence_variability(leaf_positions, orthogonal_jaw_positions,
                                                                         leaf_position_boundaries)
            result += control_point_aav * control_point_lsv * differential_weight
    return result


def small_aperture_score(beam, aperture=10):
    """ Calculate the small aperture scores of a treatment beam.

    Args:
        beam (Dataset): The beam for the small aperture score to be calculated.
        aperture (float): The aperture size to be scored.
    """
    result = 0
    previous_cumulative_weight = 0.0
    leaf_positions = []
    orthogonal_jaw_positions = []
    leaf_boundaries = _leaf_position_boundaries(beam)
    total_cumulative_weight = beam.FinalCumulativeMetersetWeight
    for control_point in beam.ControlPoints:
        differential_weight = (control_point.CumulativeMetersetWeight - previous_cumulative_weight) / \
                              total_cumulative_weight
        previous_cumulative_weight = control_point.CumulativeMetersetWeight
        control_point_leaf_positions = _control_point_leaf_positions(control_point)
        control_point_orthogonal_jaw_positions = _control_point_orthogonal_jaw_positions(control_point)
        if control_point_leaf_positions is not None:
            leaf_positions = control_point_leaf_positions
        if control_point_orthogonal_jaw_positions is not None:
            orthogonal_jaw_positions = control_point_orthogonal_jaw_positions
        if differential_weight > 0:
            result += differential_weight * _control_point_small_aperture_score(leaf_positions,
                                                                                orthogonal_jaw_positions,
                                                                                leaf_boundaries, aperture)
    return result


def _control_point_aperture_area_variability(leaf_positions, orthogonal_jaw_positions, leaf_boundaries,
                                             maximum_leaf_separations):
    leaf_pairs = int(len(leaf_positions) / 2)
    total_leaf_pair_separation = 0
    total_maximum_leaf_pair_separation = sum(maximum_leaf_separations)
    for leaf_index in range(leaf_pairs):
        if leaf_boundaries[leaf_index + 1] >= orthogonal_jaw_positions[0] \
                and leaf_boundaries[leaf_index] <= orthogonal_jaw_positions[1]:
            total_leaf_pair_separation += leaf_positions[leaf_index + leaf_pairs] - leaf_positions[leaf_index]
    if total_maximum_leaf_pair_separation > 0:
        return total_leaf_pair_separation / total_maximum_leaf_pair_separation
    return 0


def _control_point_closed_leaf_score(leaf_positions, orthogonal_jaw_positions, leaf_boundaries):
    leaf_pairs = int(len(leaf_positions) / 2)
    exposed_leaf_pairs = 0
    exposed_open_leaf_pairs = 0
    for leaf_index in range(leaf_pairs):
        if leaf_boundaries[leaf_index + 1] >= orthogonal_jaw_positions[0] \
                and leaf_boundaries[leaf_index] <= orthogonal_jaw_positions[1]:
            exposed_leaf_pairs += 1
            if leaf_positions[leaf_index + leaf_pairs] - leaf_positions[leaf_index] > 0:
                exposed_open_leaf_pairs += 1
    if exposed_leaf_pairs > 0:
        return 1 - (exposed_open_leaf_pairs / exposed_leaf_pairs)
    return 0


def _control_point_cross_axis_score(leaf_positions, orthogonal_jaw_positions, leaf_boundaries):
    leaf_pairs = int(len(leaf_positions) / 2)
    axis_crossing_leaf_pairs = 0
    exposed_open_leaf_pairs = 0
    for leaf_index in range(leaf_pairs):
        if leaf_boundaries[leaf_index + 1] >= orthogonal_jaw_positions[0] \
                and leaf_boundaries[leaf_index] <= orthogonal_jaw_positions[1]:
            if leaf_positions[leaf_index + leaf_pairs] - leaf_positions[leaf_index] > 0:
                exposed_open_leaf_pairs += 1
                if leaf_positions[leaf_index + leaf_pairs] < 0 or leaf_positions[leaf_index] > 0:
                    axis_crossing_leaf_pairs += 1
    if exposed_open_leaf_pairs > 0:
        return axis_crossing_leaf_pairs / exposed_open_leaf_pairs
    return 0


def _control_point_jaw_positions(control_point, jaw_axis='X'):
    """ Return the jaw positions for specified control point, for the corresponding jaw.

    Args:
        control_point (Dataset): The control point containing leaf position data.
    """
    for beam_limiting_device_position in control_point.BeamLimitingDevicePositions:
        if 'MLC' not in beam_limiting_device_position.RTBeamLimitingDeviceType:
            if jaw_axis in beam_limiting_device_position.RTBeamLimitingDeviceType:
                return beam_limiting_device_position.LeafJawPositions


def _control_point_leaf_positions(control_point):
    """ Return the leaf positions for specified control point.

    Args:
        control_point (Dataset): The control point containing leaf position data.
    """
    for beam_limiting_device_position in control_point.BeamLimitingDevicePositions:
        if 'MLC' in beam_limiting_device_position.RTBeamLimitingDeviceType:
            return beam_limiting_device_position.LeafJawPositions


def _control_point_leaf_sequence_variability(leaf_positions, orthogonal_jaw_positions, leaf_boundaries):
    """ Return the leaf sequence variability for control point.

     Args:
         control_point (Dataset): The control point containing leaf position data.
         leaf_boundaries (ndarray): The edge positions of leaf pairs, orthogonal to motion.
     """
    leaf_pairs = int(len(leaf_positions) / 2)
    exposed_open_leaf_pairs = 0
    leaf_positions_left = []
    leaf_positions_right = []
    delta_left = []
    delta_right = []
    for leaf_index in range(leaf_pairs - 1):
        if leaf_boundaries[leaf_index + 1] >= orthogonal_jaw_positions[0] \
                and leaf_boundaries[leaf_index] <= orthogonal_jaw_positions[1]:
            # note that paper is not clear on whether leaf pair must be open
            leaf_gap = leaf_positions[leaf_index + leaf_pairs] - leaf_positions[leaf_index]
            if leaf_gap > 0:
                exposed_open_leaf_pairs += 1
                leaf_positions_left.append(leaf_positions[leaf_index])
                leaf_positions_right.append(leaf_positions[leaf_index + leaf_pairs])
                delta_left.append(leaf_positions[leaf_index] - leaf_positions[leaf_index + 1])
                delta_right.append(leaf_positions[leaf_index + leaf_pairs] -
                                   leaf_positions[leaf_index + leaf_pairs + 1])
    min_max_left = max(leaf_positions_left) - min(leaf_positions_left)
    min_max_right = max(leaf_positions_right) - min(leaf_positions_right)
    lsv_left = sum(abs(min_max_left - np.array(delta_left)))
    lsv_right = sum(abs(min_max_right - np.array(delta_right)))
    return (lsv_left / (exposed_open_leaf_pairs * min_max_left)) * \
           (lsv_right / (exposed_open_leaf_pairs * min_max_right))


def _control_point_mean_asymmetry_distance(leaf_positions, orthogonal_jaw_positions, leaf_boundaries):
    """ Return the mean asymmetry distance for the specified control point.

     Args:
         control_point (Dataset): The control point containing leaf position data.
         leaf_boundaries (np.array): The edge positions of leaf pairs, orthogonal to motion.
     """
    leaf_pairs = int(len(leaf_positions) / 2)
    exposed_open_leaf_pairs = 0
    total_asymmetry_distance = 0
    for leaf_index in range(leaf_pairs):
        if leaf_boundaries[leaf_index + 1] >= orthogonal_jaw_positions[0] \
                and leaf_boundaries[leaf_index] <= orthogonal_jaw_positions[1]:
            leaf_gap = leaf_positions[leaf_index + leaf_pairs] - leaf_positions[leaf_index]
            leaf_sum = leaf_positions[leaf_index + leaf_pairs] + leaf_positions[leaf_index]
            if leaf_gap > 0:
                exposed_open_leaf_pairs += 1
                total_asymmetry_distance += abs(leaf_sum / 2)
    if exposed_open_leaf_pairs > 0:
        return total_asymmetry_distance / exposed_open_leaf_pairs
    return 0


def _control_point_orthogonal_jaw_positions(control_point):
    """ Return the jaw positions for specified control point, for the jaws orthogonal to leaves.

    Args:
        control_point (Dataset): The control point containing leaf position data.
    """
    orthogonal = {'X': 'Y', 'Y': 'X'}
    for beam_limiting_device_position in control_point.BeamLimitingDevicePositions:
        if 'MLC' in beam_limiting_device_position.RTBeamLimitingDeviceType:
            mlc_type = beam_limiting_device_position.RTBeamLimitingDeviceType.replace('MLC', '')
            return _control_point_jaw_positions(control_point, orthogonal[mlc_type])


def _control_point_running_jaw_positions(control_point):
    """ Return the jaw positions for specified control point, for the jaws parallel to leaves.

    Args:
        control_point (Dataset): The control point containing leaf position data.
    """
    running = {'X': 'X', 'Y': 'Y'}
    for beam_limiting_device_position in control_point.BeamLimitingDevicePositions:
        if 'MLC' in beam_limiting_device_position.RTBeamLimitingDeviceType:
            mlc_type = beam_limiting_device_position.RTBeamLimitingDeviceType.replace('MLC', '')
            return _control_point_jaw_positions(control_point, running[mlc_type])


def _control_point_small_aperture_score(leaf_positions, orthogonal_jaw_positions, leaf_boundaries, aperture):
    """ Return the number of open leaf pair separations less than specified aperture threshold.

    Args:
        leaf_positions (ndarray): The leaf positions for the control point.
        orthogonal_jaw_positions (ndarray): The positions of the orthogonal jaws.
        leaf_boundaries (ndarray): The edge positions of leaf pairs, orthogonal to motion.
        aperture (float): The threshold defining a small aperture.
    """
    leaf_pairs = int(len(leaf_positions) / 2)
    exposed_open_leaf_pairs = 0
    small_aperture_leaf_pairs = 0
    for leaf_index in range(leaf_pairs):
        if leaf_boundaries[leaf_index + 1] >= orthogonal_jaw_positions[0] \
                and leaf_boundaries[leaf_index] <= orthogonal_jaw_positions[1]:
            leaf_gap = leaf_positions[leaf_index + leaf_pairs] - leaf_positions[leaf_index]
            if leaf_gap > 0:
                exposed_open_leaf_pairs += 1
                if leaf_gap < aperture:
                    small_aperture_leaf_pairs += 1
    if exposed_open_leaf_pairs > 0:
        return small_aperture_leaf_pairs / exposed_open_leaf_pairs
    return 0


def _control_point_weight(beam, control_point_index):
    """ Determine the fractional differential meterset weight for specified control point.

    Args:
        beam (Dataset): The beam containing the control point for which weight is to be calculated.
        control_point_index (int): The index of the control point dataset in the control point sequence.
    """
    if control_point_index == 0:
        return beam.ControlPoints[0].CumulativeMetersetWeight / beam.FinalCumulativeMetersetWeight
    else:
        return (beam.ControlPoints[control_point_index].CumulativeMetersetWeight -
                beam.ControlPoints[control_point_index - 1].CumulativeMetersetWeight) \
               / beam.FinalCumulativeMetersetWeight


def _has_multi_leaf_collimator(beam):
    """ Determines if beam has multi-leaf collimator device type.

    Args:
        beam (Dataset): The beam possibly containing multi-leaf collimator device type.
    """
    for beam_limiting_device in beam.BeamLimitingDevices:
        if 'MLC' in beam_limiting_device.RTBeamLimitingDeviceType:
            return True
    return False


def _leaf_position_boundaries(beam):
    """ Determine leaf position boundaries for MLC system specified in plan

    Args:
        beam (Dataset): The beam for the small aperture score to be calculated.
    """
    for beam_limiting_device in beam.BeamLimitingDevices:
        if 'MLC' in beam_limiting_device.RTBeamLimitingDeviceType:
            return beam_limiting_device.LeafPositionBoundaries


def _maximum_leaf_separations(beam):
    """ Return the maximum leaf separations for each leaf pair over the entire beam.

    Args:
        beam (Dataset): The beam containing control points with leaf positions.
    """
    # pre-conditions
    if not isinstance(beam, dicom.dataset.Dataset):
        raise TypeError('Beam is not valid type.')
    if not _has_multi_leaf_collimator(beam):
        raise ValueError('Beam has no multi-leaf collimator.')
    number_leaf_pairs = len(_leaf_position_boundaries(beam)) - 1
    leaf_separations = np.array([0] * number_leaf_pairs)
    leaf_pair_indices = range(number_leaf_pairs)
    for control_point in beam.ControlPoints:
        leaf_positions = _control_point_leaf_positions(control_point)
        for index in leaf_pair_indices:
            curr_aperture = leaf_positions[index + number_leaf_pairs] - leaf_positions[index]
            if curr_aperture > leaf_separations[index]:
                leaf_separations[index] = curr_aperture
    return leaf_separations


def _weighted_metric(beam, metric):
    """ Return the MU weighted average of a control point metric.

    Args:
        beam (Dataset): The beam containing control points.
        metric (string): The metric to be calculated.
    """
    result = 0
    previous_cumulative_weight = 0.0
    total_cumulative_weight = beam.FinalCumulativeMetersetWeight
    leaf_positions = []
    orthogonal_jaw_positions = []
    leaf_boundaries = _leaf_position_boundaries(beam)
    for control_point in beam.ControlPoints:
        differential_weight = (control_point.CumulativeMetersetWeight - previous_cumulative_weight) / \
                              total_cumulative_weight
        previous_cumulative_weight = control_point.CumulativeMetersetWeight
        control_point_leaf_positions = _control_point_leaf_positions(control_point)
        control_point_orthogonal_jaw_positions = _control_point_orthogonal_jaw_positions(control_point)
        if control_point_leaf_positions is not None:
            leaf_positions = control_point_leaf_positions
        if control_point_orthogonal_jaw_positions is not None:
            orthogonal_jaw_positions = control_point_orthogonal_jaw_positions
        if differential_weight > 0:
            result += differential_weight * _control_point_metric_list[metric](leaf_positions,
                                                                               orthogonal_jaw_positions,
                                                                               leaf_boundaries)
    return result

# define metric dictionary
_metric_list = {'AAV': aperture_area_variability,
                'CLS': closed_leaf_score,
                'CAS': cross_axis_score,
                'FMC': fluence_map_complexity,
                'LSV': leaf_sequence_variability,
                'MAD': mean_asymmetry_distance,
                'MCS': modulation_complexity_score}

# define control point metric dictionary
_control_point_metric_list = {'CLS': _control_point_closed_leaf_score,
                              'CAS': _control_point_cross_axis_score,
                              'LSV': _control_point_leaf_sequence_variability,
                              'MAD': _control_point_mean_asymmetry_distance}