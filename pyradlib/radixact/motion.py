# -*- coding: utf-8 -*-
""" Motion analysis module.

This module provides functionality for processing of motion data from Synchrony treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import datetime
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import os
import pandas as pd
import re
import xml.etree.ElementTree as et


def read_motion_data(xml_filepaths: npt.ArrayLike):
    """Extracts motion data from *motionData.xml files representing a single fraction,
    as produced for the Radixact Delivery Analysis tool.

    Parameters
    ----------
    xml_filepaths : array_like
        List of paths to XML files containing motion data.

    Returns
    -------
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset : array_like
        Output arrays, containing paired timestamps, potential differences, rigid body
        differences, and IEC-X, IEC-Y and IEC-Z offsets (in sec, mm).

    Notes
    -----
    These files are cached in C:/tomo/da/pts/URnumber/*motionData.xml when patient data
    is loaded within the Delivery Analysis tool. There may be multiple files per fraction.

    Pauses in fractions (e.g. due to user, Synchrony, or system errors) are indicated by
    numpy.nan. This allows more accurate plotting of the data using matplotlib.
    """
    time = []
    potential_diff = []
    rigid_body = []
    x_offset = []
    y_offset = []
    z_offset = []
    for xml_path in xml_filepaths:
        tree = et.parse(xml_path)
        root = tree.getroot()
        for delivery_segment in root.iter("DeliverySegment"):
            for radiation_results in delivery_segment.iter("RadiationResults"):
                for model_data in radiation_results.iter("ModelData"):
                    for datapoint in model_data.iter("DataPoint"):
                        time.append(float(datapoint[0].text) / 1e3)
                        potential_diff.append(float(datapoint[1].text))
                        rigid_body.append(float(datapoint[2].text))
                        x_offset.append(float(datapoint[3][0].text))
                        y_offset.append(float(datapoint[3][1].text))
                        z_offset.append(float(datapoint[3][2].text))
                    # include NaN to reflect beam being paused at end of radiation result (which forces matplotlib to break line)
                    time.append(np.nan)
                    potential_diff.append(np.nan)
                    rigid_body.append(np.nan)
                    x_offset.append(np.nan)
                    y_offset.append(np.nan)
                    z_offset.append(np.nan)
    # remove trailing NaN, because completion is not a pause
    time.pop()
    potential_diff.pop()
    rigid_body.pop()
    x_offset.pop()
    y_offset.pop()
    z_offset.pop()
    return time, potential_diff, rigid_body, x_offset, y_offset, z_offset


def convert_motion_data_to_csv(xml_filepaths: npt.ArrayLike, csv_filepath: str):
    """Converts motion data from *motionData.xml files produced for the Radixact
    Delivery Analysis tool to CSV files.

    Parameters
    ----------
    xml_filepaths : array_like
        List of paths to XML files containing motion data.
    csv_filepath : str
        Path for CSV file to be written.
    """
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset = read_motion_data(
        xml_filepaths
    )
    np.savetxt(
        csv_filepath,
        zip(time, x_offset, y_offset, z_offset),
        delimiter=",",
        header="Time (s), IEC-X (mm), IEC-Y (mm), IEC-Z (mm)",
    )


def modify_motion_data(
    time: npt.ArrayLike,
    potential_diff: npt.ArrayLike,
    rigid_body: npt.ArrayLike,
    x_offset: npt.ArrayLike,
    y_offset: npt.ArrayLike,
    z_offset: npt.ArrayLike,
    remove_duplicates: bool = True,
    remove_variable_pauses: bool = True,
    replacement_pause_length: float = 5,
    zero_reference_point: bool = False,
):
    """Modify extracted motion data for improved presentation.

    Parameters
    ----------
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset : array_like
        Input arrays, containing paired timestamps, potential differences, rigid body
        differences, and IEC-X, IEC-Y and IEC-Z offsets (in sec, mm).
    remove_duplicates : bool, optional
        Flag indicating removal of duplicate data points, including at start of
        treatment, for reduced data density and profile length (default is True).
    remove_variable_pauses : bool, optional
        Flag indicating replacement of variable length pauses in treatment with short
        fixed length, defined by pause_offset_time (default is True).
    replacement_pause_length : float, optional
        Time period used when replacing variable pauses in beam delivery, e.g.,
        5 sec (default).
    zero_reference_point : bool, optional
        Flag indicating whether coordinates are shifted for starting positions, i.e.
        whether X, Y and Z are shifted to 0 at time 0, and following any pause (default
        is False).

    Returns
    -------
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset : array_like
        Output arrays, containing adjusted paired timestamps, potential differences,
        rigid body differences, and IEC-X, IEC-Y and IEC-Z offsets (in sec, mm).
    pause_times : array_like
        Output array, containing 0 or more paired timestamps indicating when pauses in
        treatment occur, for plotting (in sec, mm).
    """
    # adjust time and position data
    deleted_indices = []
    if remove_duplicates:
        for i in range(len(time) - 1):
            if (
                np.array([x_offset[i], y_offset[i], z_offset[i]])
                == np.array([x_offset[i + 1], y_offset[i + 1], z_offset[i + 1]])
            ).all():
                deleted_indices.append(i)
    new_time = np.delete(time, deleted_indices)
    new_potential_diff = np.delete(potential_diff, deleted_indices)
    new_rigid_body = np.delete(rigid_body, deleted_indices)
    new_x_offset = np.delete(x_offset, deleted_indices)
    new_y_offset = np.delete(y_offset, deleted_indices)
    new_z_offset = np.delete(z_offset, deleted_indices)
    # adjust time data
    adjusted_time = []
    pauses = []
    # if no pauses, shift so first image is time 0
    if np.isnan(time).sum() == 0:
        # if no pauses, shift so first result is time 0
        adjusted_time = np.array(new_time) - new_time[0]
    else:
        # if there are pauses, contract and/or calculate pause lengths
        for i in range(len(new_time)):
            if i == 0:
                adjusted_time.append(0)
            elif np.isnan(time[i]):
                adjusted_time.append(np.nan)
            else:
                if np.isnan(adjusted_time[-1]) and remove_variable_pauses:
                    last_valid_index = np.where(~np.isnan(adjusted_time))[0][-1]
                    adjusted_time.append(
                        adjusted_time[last_valid_index] + replacement_pause_length
                    )
                    pauses.append(
                        [
                            adjusted_time[last_valid_index],
                            adjusted_time[last_valid_index] + replacement_pause_length,
                        ]
                    )
                elif np.isnan(adjusted_time[-1]) and not remove_variable_pauses:
                    last_valid_index = np.where(~np.isnan(adjusted_time))[0][-1]
                    adjusted_time.append(
                        adjusted_time[last_valid_index]
                        + (new_time[i] - new_time[last_valid_index])
                    )
                    pauses.append([adjusted_time[last_valid_index], adjusted_time[-1]])
                else:
                    adjusted_time.append(
                        adjusted_time[-1] + (new_time[i] - new_time[i - 1])
                    )
    if zero_reference_point:
        x_shift = np.NaN
        y_shift = np.NaN
        z_shift = np.NaN
        for i in range(len(adjusted_time)):
            if np.isnan(new_x_offset[i]):
                # if pause encountered
                x_shift = np.NaN
                y_shift = np.NaN
                z_shift = np.NaN
            else:
                if np.isnan(x_shift):
                    # if at start of treatment or following pause
                    x_shift = -new_x_offset[i]
                    y_shift = -new_y_offset[i]
                    z_shift = -new_z_offset[i]
                new_x_offset[i] = new_x_offset[i] + x_shift
                new_y_offset[i] = new_y_offset[i] + y_shift
                new_z_offset[i] = new_z_offset[i] + z_shift
    return (
        adjusted_time,
        new_potential_diff,
        new_rigid_body,
        new_x_offset,
        new_y_offset,
        new_z_offset,
        pauses,
    )


def patient_fraction_xml_lists(patient_path: str):
    """Produces list containing lists of *motionData.xml paths for each fraction of one or
    more treatments for a given patient.

    Parameters
    ----------
    patient_path : str
        The path containing the *motionData.xml files to be listed.

    Notes
    -----
    The specified directory should correspond with those cached in C:/tomo/da/pts/URnumber/
    when patient data is loaded within the Delivery Analysis tool.
    """
    motion_files = sorted(
        glob.glob(os.path.join(patient_path, "*motionData.xml")),
        key=lambda x: float(
            re.findall("(\d+.\d+)", x.split("-")[-2] + "." + x.split("-")[-1])[0]
        ),
    )
    if len(motion_files) == 0:
        return []
    motion_paths = {}
    for motion_file in motion_files:
        id = (
            "-".join(os.path.split(motion_file)[1].split("-")[0:-2])
            + "-"
            + os.path.split(motion_file)[1].split("-")[-2].zfill(2)
        )
        if id in motion_paths:
            curr_list = motion_paths[id]
            curr_list.append(motion_file)
            motion_paths[id] = curr_list
        else:
            motion_paths[id] = [motion_file]
    motion_paths = dict(sorted(motion_paths.items()))
    return list(motion_paths.values())


def convert_patient_motion_data_to_csv(patient_path: str, csv_path: str):
    fraction_xml_paths = patient_fraction_xml_lists(patient_path)
    for fraction in range(len(fraction_xml_paths)):
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset = (
            read_motion_data(fraction_xml_paths[fraction])
        )
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset, pauses = (
            modify_motion_data(
                time,
                potential_diff,
                rigid_body,
                x_offset,
                y_offset,
                z_offset,
                remove_duplicates=True,
                remove_variable_pauses=False,
                # replacement_pause_length =
                zero_reference_point=True,
            )
        )
        np.savetxt(
            os.path.join(csv_path, "Fraction " + str(fraction + 1) + ".csv"),
            np.column_stack((time, x_offset, y_offset, z_offset)),
            delimiter=",",
            header="Time (s), IEC-X (mm), IEC-Y (mm), IEC-Z (mm)",
        )


def calculate_vector_displacements(
    time: npt.ArrayLike,
    x_offset: npt.ArrayLike,
    y_offset: npt.ArrayLike,
    z_offset: npt.ArrayLike,
    point_by_point_displacement: bool = False,
):
    """Calculate 3D vector displacements.

    Parameters
    ----------
    time, x_offset, y_offset, z_offset : array_like
        Input arrays, containing paired timestamps, and IEC-X, IEC-Y and IEC-Z offsets
        (in sec, mm, with duplicates removed).
    point_by_point_displacement : bool, optional
        Flag indicating whether displacement should be calculated relative to most recent
        offset data, as opposed to most recent starting offset data (i.e., at start of
        treatment or after a pause, default is False).

    Returns
    -------
    vector_displacement_times, vector_displacements : array_like
        Output arrays, containing adjusted paired timestamps, and displacement times
        (in sec, mm).
    """
    result_time = []
    result_displacement = []
    x_ref = np.NaN
    y_ref = np.NaN
    z_ref = np.NaN
    for i in range(len(time)):
        if np.isnan(x_offset[i]):
            # if pause encountered
            x_ref = np.NaN
            y_ref = np.NaN
            z_ref = np.NaN
        else:
            if np.isnan(x_ref):
                # if reference not set (i.e., at start, or following pause)
                x_ref = x_offset[i]
                y_ref = x_offset[i]
                z_ref = x_offset[i]
            else:
                result_time.append(time[i])
                result_displacement.append(
                    np.sqrt(
                        np.square(x_ref - x_offset[i])
                        + np.square(y_ref - y_offset[i])
                        + np.square(z_ref - z_offset[i])
                    )
                )
                if point_by_point_displacement:
                    x_ref = x_offset[i]
                    y_ref = x_offset[i]
                    z_ref = x_offset[i]
    return np.array(result_time), np.array(result_displacement)


def plot_motion_data(
    fraction_xml_paths: npt.ArrayLike,
    png_filepath: str,
    title: str,
    plot_potential_diff: bool = False,
    plot_rigid_body: bool = False,
    number_columns: int = 5,
    share_x_axis: bool = True,
    remove_duplicates: bool = True,
    remove_variable_pauses: bool = True,
    replacement_pause_length: float = 5,
    pause_color="lightgrey",
    zero_reference_point: bool = False,
):
    """Plot collection of motion data spanning multiple fractions.

    Parameters
    ----------
    fraction_xml_paths : array_like
        List containing lists of *motionData.xml paths for each fraction.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.
    plot_potential_diff : bool, optional
        Flag indicating whether potential difference is plotted (default is False).
    plot_rigid_body : bool, optional
        Flag indicating whether rigid body difference is plotted (default is False).
    number_columns : int, optional
        Number of columns for plotted data.
    remove_duplicates : bool, optional
        Flag indicating removal of duplicate data points, including at start of treatment,
        for reduced data density and profile length (default is True).
    remove_variable_pauses : bool, optional
        Flag indicating replacement of variable length pauses in treatment with short
        fixed length, defined by pause_offset_time (default is True).
    replacement_pause_length : float, optional
        Time period used when replacing variable pauses in beam delivery, e.g., 5 sec
        (default).
    pause_colour : str, optional
        Color to use to indicate pauses in treatment, e.g., 'lightgrey' (default).
    zero_reference_point : bool, optional
        Flag indicating whether coordinates are shifted for starting positions, i.e.
        whether X, Y and Z are shifted to 0 at time 0, and following any pause (default
        is False).
    """
    ncols = np.min([number_columns, len(fraction_xml_paths)])
    nrows = int((len(fraction_xml_paths) + (ncols - 1)) / ncols)
    fig, axs = plt.subplots(
        ncols=ncols,
        nrows=nrows,
        sharex=share_x_axis,
        sharey=True,
        gridspec_kw={"hspace": 0, "wspace": 0},
        constrained_layout=True,
        figsize=(ncols * 2 + 1, nrows * 2 + 1),
    )
    for fraction in range(len(fraction_xml_paths)):
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset = (
            read_motion_data(fraction_xml_paths[fraction])
        )
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset, pauses = (
            modify_motion_data(
                time,
                potential_diff,
                rigid_body,
                x_offset,
                y_offset,
                z_offset,
                remove_duplicates,
                remove_variable_pauses,
                replacement_pause_length,
                zero_reference_point,
            )
        )
        col = int(fraction % ncols)
        row = int(fraction / ncols)
        if nrows > 1:
            (l1,) = axs[row, col].plot(time, x_offset)
            (l2,) = axs[row, col].plot(time, y_offset)
            (l3,) = axs[row, col].plot(time, z_offset)
        elif ncols > 1:
            (l1,) = axs[col].plot(time, x_offset)
            (l2,) = axs[col].plot(time, y_offset)
            (l3,) = axs[col].plot(time, z_offset)
        else:
            (l1,) = axs.plot(time, x_offset)
            (l2,) = axs.plot(time, y_offset)
            (l3,) = axs.plot(time, z_offset)
        plotted = [l1, l2, l3]
        legend = ["IEC-X", "IEC-Y", "IEC-Z"]
        if plot_potential_diff:
            if nrows > 1:
                (l4,) = axs[row, col].plot(time, potential_diff)
            elif ncols > 1:
                (l4,) = axs[col].plot(time, potential_diff)
            else:
                (l4,) = axs.plot(time, potential_diff)
            plotted.append(l4)
            legend.append("Potential Diff")
        if plot_rigid_body:
            if nrows > 1:
                (l5,) = axs[row, col].plot(time, rigid_body)
            elif ncols > 1:
                (l5,) = axs[col].plot(time, rigid_body)
            else:
                (l5,) = axs.plot(time, rigid_body)
            plotted.append(l5)
            legend.append("Rigid Body")
        for pause_time in pauses:
            if nrows > 1:
                axs[row, col].axvspan(
                    pause_time[0], pause_time[1], color=pause_color, lw=0
                )
            elif ncols > 1:
                axs[col].axvspan(pause_time[0], pause_time[1], color=pause_color, lw=0)
            else:
                axs.axvspan(pause_time[0], pause_time[1], color=pause_color, lw=0)
        if nrows > 1:
            axs[row, col].set_title(fraction + 1, y=1.0, pad=-14)
        elif ncols > 1:
            axs[col].set_title(fraction + 1, y=1.0, pad=-14)
        else:
            axs.set_title(fraction + 1, y=1.0, pad=-14)
    # plt.ylim(-_plot_maximum_motion,_plot_maximum_motion)
    fig.legend(
        plotted,
        legend,
        loc="center right",
        ncol=ncols,
        fancybox=True,
        bbox_to_anchor=(1.05, 1),
    )
    fig.supxlabel("Time (s)")
    fig.supylabel("Displacement (mm)")
    # remove empty grids
    if len(fraction_xml_paths) % ncols > 0:
        if nrows > 1:
            axs[-1, -1].set_axis_off()
            axs[-2, -1].xaxis.set_tick_params(which="both", labelbottom=True)
        elif ncols > 1:
            axs[-1].set_axis_off()
        if len(fraction_xml_paths) % ncols < 4:
            if nrows > 1:
                axs[-1, -2].set_axis_off()
                axs[-2, -2].xaxis.set_tick_params(which="both", labelbottom=True)
            elif ncols > 1:
                axs[-2].set_axis_off()
            if len(fraction_xml_paths) % ncols < 3:
                if nrows > 1:
                    axs[-1, -3].set_axis_off()
                    axs[-2, -3].xaxis.set_tick_params(which="both", labelbottom=True)
                elif ncols > 1:
                    axs[-3].set_axis_off()
                if len(fraction_xml_paths) % ncols < 2:
                    if nrows > 1:
                        axs[-1, -4].set_axis_off()
                        axs[-2, -4].xaxis.set_tick_params(
                            which="both", labelbottom=True
                        )
                    elif ncols > 1:
                        axs[-4].set_axis_off()
    fig.suptitle(title)
    plt.savefig(png_filepath, bbox_inches="tight")


def plot_patient_motion_data(
    patient_path: str, png_filepath: str, title: str, zero_reference_point: bool = False
):
    """Plot motion data contained within specific patient directory.

    Parameters
    ----------
    patient_path : str
        The path containing the *motionData.xml files to be plotted.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.
    zero_reference_point : bool, optional
        Flag indicating whether coordinates are shifted for starting positions, i.e.
        whether X, Y and Z are shifted to 0 at time 0, and following any pause (default
        is False).

    Notes
    -----
    The specified directory should correspond with those cached in C:/tomo/da/pts/URnumber/
    when patient data is loaded within the Delivery Analysis tool.
    """
    motion_paths = patient_fraction_xml_lists(patient_path)
    if len(motion_paths) > 0:
        plot_motion_data(
            motion_paths,
            png_filepath,
            title,
            zero_reference_point=zero_reference_point,
        )


def plot_patient_motion_data_boxplot(
    patient_path: str,
    png_filepath: str,
    title: str,
    zero_reference_point: bool = True,
    type="classical",
):
    """Plot boxplot of motion data contained within specific patient directory.

    Parameters
    ----------
    patient_path : str
        The path containing the *motionData.xml files to be plotted.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.
    zero_reference_point : bool, optional
        Flag indicating whether coordinates are shifted for starting positions, i.e.
        whether X, Y and Z are shifted to 0 at time 0, and following any pause (default
        is False).
    type : str, optional
        Type of boxplot to produce, possible values are "classical" and "functional" 
        (default is "classical").
    """
    motion_paths = patient_fraction_xml_lists(patient_path)
    if len(motion_paths) == 0:
        return
    fig, axs = plt.subplots(
        ncols=1,
        nrows=3,
        sharex=True,
        sharey=True,
        gridspec_kw={"hspace": 0, "wspace": 0},
        constrained_layout=True,
        figsize=(4, 7),
    )
    all_x_offset = []
    all_y_offset = []
    all_z_offset = []
    for fraction in range(len(motion_paths)):
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset = (
            read_motion_data(motion_paths[fraction])
        )
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset, pauses = (
            modify_motion_data(
                time,
                potential_diff,
                rigid_body,
                x_offset,
                y_offset,
                z_offset,
                remove_duplicates=True,
                remove_variable_pauses=True,
                replacement_pause_length=0,
                zero_reference_point=True,
            )
        )
        all_x_offset.append(x_offset[~np.isnan(x_offset)])
        all_y_offset.append(y_offset[~np.isnan(y_offset)])
        all_z_offset.append(z_offset[~np.isnan(z_offset)])
    if type == "classical":
        axs[0].boxplot(all_x_offset)
        axs[1].boxplot(all_y_offset)
        axs[2].boxplot(all_z_offset)
    if type == "functional":
        plots = []
        for level in [100, 50]:
            color = mpl.colormaps["Blues"](int(25.6 + ((100 - level) * 2)))
            plots.append(
                axs[0].fill_between(
                    range(len(all_x_offset)),
                    [np.percentile(data, level) for data in all_x_offset],
                    [np.percentile(data, 100 - level) for data in all_x_offset],
                    color=color,
                )
            )
            axs[1].fill_between(
                range(len(all_y_offset)),
                [np.percentile(data, level) for data in all_y_offset],
                [np.percentile(data, 100 - level) for data in all_y_offset],
                color=color,
            )
            axs[2].fill_between(
                range(len(all_z_offset)),
                [np.percentile(data, level) for data in all_z_offset],
                [np.percentile(data, 100 - level) for data in all_z_offset],
                color=color,
            )
        (l1,) = axs[0].plot(
            range(len(all_x_offset)),
            [np.median(data) for data in all_x_offset],
            color=mpl.colormaps["Blues"](226),
        )
        plots.append(l1)
        axs[1].plot(
            range(len(all_x_offset)),
            [np.median(data) for data in all_y_offset],
            color=mpl.colormaps["Blues"](226),
        )
        axs[2].plot(
            range(len(all_x_offset)),
            [np.median(data) for data in all_z_offset],
            color=mpl.colormaps["Blues"](226),
        )
    axs[0].set_title("IEC-X", y=1.0, pad=-14)
    axs[1].set_title("IEC-Y", y=1.0, pad=-14)
    axs[2].set_title("IEC-Z", y=1.0, pad=-14)
    axs[2].tick_params(labelbottom=False)
    fig.supxlabel("Fraction")
    fig.supylabel("Displacement (mm)")
    fig.suptitle(title)
    plt.savefig(png_filepath, bbox_inches="tight")


def has_motion_data(patient_path: str):
    """Check whether motion data exists within specific patient directory.

    Parameters
    ----------
    patient_path : str
        The path nominally containing *motionData.xml files.

    Returns
    -------
    result : bool
        True, if *motionData.xml found, otherwise False.
    """
    files = glob.glob(os.path.join(patient_path, "*motionData.xml"))
    return len(files) > 0


def plot_vector_displacement_probability(
    fraction_xml_paths: npt.ArrayLike,
    png_filepath: str,
    title: str,
    offset_max: float = 5,
    offset_bin: float = 0.25,
    vector_max: float = 10,
    vector_bin: float = 0.25,
    stacked_colour_histogram: bool = False,
):
    """Plot vector displacement probability of motion data spanning multiple fractions.

    Parameters
    ----------
    fraction_xml_paths : array_like
        List containing lists of *motionData.xml paths for each fraction of one or
        more treatments.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.
    offset_max, vector_max : float
        Displacement limit for plotting for IEC axes and vector (in mm, 5 and 10 default).
    offset_bin, vector_bin : float
        Width of displacement bin for histograms for IEC axes and vector (in mm, 0.25 default).
    stacked_colour_histogram : bool
        Flag to indicate whether to stack histogram bars according to time of displacement.
    """
    time_combined = []
    x_displacements = []
    y_displacements = []
    z_displacements = []
    disp_time_combined = []
    v_displacements = []
    for fraction in range(len(fraction_xml_paths)):
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset = (
            read_motion_data(fraction_xml_paths[fraction])
        )
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset, pauses = (
            modify_motion_data(
                time,
                potential_diff,
                rigid_body,
                x_offset,
                y_offset,
                z_offset,
                zero_reference_point=True,
            )
        )
        disp_time, disp = calculate_vector_displacements(
            time, x_offset, y_offset, z_offset
        )
        time_combined = [*time_combined, *time]
        disp_time_combined = [*disp_time_combined, *disp_time]
        x_displacements = [*x_displacements, *x_offset]
        y_displacements = [*y_displacements, *y_offset]
        z_displacements = [*z_displacements, *z_offset]
        v_displacements = [*v_displacements, *disp]
    time_combined = np.array(time_combined)
    disp_time_combined = np.array(disp_time_combined)
    x_displacements = np.array(x_displacements)
    y_displacements = np.array(y_displacements)
    z_displacements = np.array(z_displacements)
    v_displacements = np.array(v_displacements)
    label = ["Displacement"]
    if stacked_colour_histogram:
        max_time = np.nanmax(time_combined)
        num_stacks = int(max_time / 100) + 1
        x_split = []
        y_split = []
        z_split = []
        v_split = []
        label = []
        for stack in range(num_stacks):
            mask_disp = (time_combined >= (stack * 100)) & (
                time_combined < (stack * 100 + 100)
            )
            mask_vect = (disp_time_combined >= (stack * 100)) & (
                disp_time_combined < (stack * 100 + 100)
            )
            x_split.append(x_displacements[mask_disp])
            y_split.append(y_displacements[mask_disp])
            z_split.append(z_displacements[mask_disp])
            v_split.append(v_displacements[mask_vect])
            label.append(str(stack * 100) + "-" + str(stack * 100 + 100) + " s")
        x_displacements = x_split
        y_displacements = y_split
        z_displacements = z_split
        v_displacements = v_split
    fig, axs = plt.subplots(
        ncols=4,
        nrows=1,
        gridspec_kw={"hspace": 0, "wspace": 0},
        constrained_layout=True,
        figsize=(9, 3),
    )
    axs[0].hist(
        x_displacements,
        list(np.arange(-offset_max, offset_max + offset_bin, offset_bin)),
        density=False,
        stacked=stacked_colour_histogram,
    )
    if not stacked_colour_histogram:
        axs[0].axvline(np.nanmean(x_displacements), color="r")
    axs[0].set_title("IEC-X")
    axs[1].hist(
        y_displacements,
        list(np.arange(-offset_max, offset_max + offset_bin, offset_bin)),
        density=False,
        stacked=stacked_colour_histogram,
    )
    if not stacked_colour_histogram:
        axs[1].axvline(np.nanmean(y_displacements), color="r")
    axs[1].set_title("IEC-Y")
    axs[2].hist(
        z_displacements,
        list(np.arange(-offset_max, offset_max + offset_bin, offset_bin)),
        density=False,
        stacked=stacked_colour_histogram,
    )
    if not stacked_colour_histogram:
        axs[2].axvline(np.nanmean(z_displacements), color="r")
    axs[2].set_title("IEC-Z")
    axs[3].hist(
        v_displacements,
        list(np.arange(0, vector_max + vector_bin, vector_bin)),
        density=False,
        cumulative=-1,
        stacked=stacked_colour_histogram,
    )
    axs[3].set_title("Vector")
    if not stacked_colour_histogram:
        axs[3].axvline(np.percentile(v_displacements, 95), color="r")
        axs[3].text(
            np.percentile(v_displacements, 95) + 0.5,
            0.75 * len(v_displacements),
            "$r_{95}$ = \n"
            + "{:0.1f}".format(np.percentile(v_displacements, 95))
            + "mm",
            color="r",
        )
    fig.supylabel("Number of data points")
    fig.supxlabel("Displacement (mm)")
    fig.suptitle(title)
    if stacked_colour_histogram:
        fig.legend(
            label,
            title="Time period",
            loc="center right",
            fancybox=True,
            bbox_to_anchor=(1.15, 0.5),
        )
    plt.savefig(png_filepath, bbox_inches="tight")


def plot_patient_vector_displacement_probability(
    patient_path: str, png_path: str, title: str, stacked_colour_histogram: bool = False
):
    """Plot vector displacement probability for motion data contained within specific
    patient directory.

    Parameters
    ----------
    patient_path : str
        The path containing the *motionData.xml files to be plotted.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.
    stacked_colour_histogram : bool
        Flag to indicate whether to stack histogram bars according to time of displacement.

    Notes
    -----
    The specified directory should correspond with those cached in C:/tomo/da/pts/URnumber/
    when patient data is loaded within the Delivery Analysis tool.
    """
    motion_paths = patient_fraction_xml_lists(patient_path)
    plot_vector_displacement_probability(
        motion_paths,
        png_path,
        title,
        stacked_colour_histogram=stacked_colour_histogram,
    )


def plot_cohort_vector_displacement_probability(
    cohort_path: str, png_path: str, title: str, stacked_colour_histogram: bool = False
):
    """Plot vector displacement probability for motion data contained within cohort
    directory, where motion data is nested in patient and Delivery Analysis folders.

    Parameters
    ----------
    cohort_path : str
        The path containing patient directories, with *motionData.xml files in XML specific
        sub-directories, to be plotted.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.
    stacked_colour_histogram : bool
        Flag to indicate whether to stack histogram bars according to time of displacement.

    Notes
    -----
    The specified directory should correspond two levels higher than those cached in
    C:/tomo/da/pts/URnumber/ when patient data is loaded within the Delivery Analysis tool.
    """
    motion_files = sorted(
        glob.glob(os.path.join(cohort_path, "**/**/*motionData.xml")),
        key=lambda x: float(
            re.findall("(\d+.\d+)", x.split("-")[-2] + "." + x.split("-")[-1])[0]
        ),
    )
    motion_paths = {}
    for motion_file in motion_files:
        id = (
            "-".join(os.path.split(motion_file)[1].split("-")[0:-2])
            + "-"
            + os.path.split(motion_file)[1].split("-")[-2].zfill(2)
        )
        if id in motion_paths:
            curr_list = motion_paths[id]
            curr_list.append(motion_file)
            motion_paths[id] = curr_list
        else:
            motion_paths[id] = [motion_file]
    motion_paths = dict(sorted(motion_paths.items()))
    plot_vector_displacement_probability(
        list(motion_paths.values()),
        png_path,
        title,
        stacked_colour_histogram=stacked_colour_histogram,
    )


def analyse_motion_data(fraction_xml_paths: npt.ArrayLike, all_fraction_label="all"):
    """Analyse motion data spanning multiple fractions.

    Parameters
    ----------
    fraction_xml_paths : array_like
        List containing lists of *motionData.xml paths for each fraction of one or
        more treatments.
    all_fraction_label : str
        Text to be used in fraction column to indicate summary of treatment results.

    Returns
    -------
    results : Pandas DataFrame
        A dataframe containing data extracted from the XML files.
    """
    results = []
    time_combined = []
    potential_diff_combined = []
    rigid_body_combined = []
    x_displacements = []
    y_displacements = []
    z_displacements = []
    disp_time_combined = []
    v_displacements = []
    header = [
        "Fraction",
        "Data acquisition duration",
        "Number of delivery fragments",
        "Pause duration",
        "Active duration",
        "IEC-X mean",
        "IEC-X standard deviation",
        "IEC-X median",
        "IEC-X range",
        "IEC-Y mean",
        "IEC-Y standard deviation",
        "IEC-Y median",
        "IEC-Y range",
        "IEC-Z mean",
        "IEC-Z standard deviation",
        "IEC-Z median",
        "IEC-Z range",
        "Vector displacement mean",
        "Vector displacement standard deviation",
        "Vector displacement median",
        "Vector displacement 80th percentile",
        "Vector displacement 90th percentile",
        "Vector displacement 95th percentile",
        "Vector displacement maximum",
        "Rigid body mean",
        "Rigid body standard deviation",
        "Rigid body median",
        "Rigid body maximum",
    ]

    def extract_metrics(x_data, y_data, z_data, disp_data, rigid_body_data):
        metrics = []
        metrics.append(np.nanmean(x_data))
        metrics.append(np.nanstd(x_data))
        metrics.append(np.nanmedian(x_data))
        metrics.append(np.nanmax(x_data) - np.nanmin(x_data))
        metrics.append(np.nanmean(y_data))
        metrics.append(np.nanstd(y_data))
        metrics.append(np.nanmedian(y_data))
        metrics.append(np.nanmax(y_data) - np.nanmin(y_data))
        metrics.append(np.nanmean(z_data))
        metrics.append(np.nanstd(z_data))
        metrics.append(np.nanmedian(z_data))
        metrics.append(np.nanmax(z_data) - np.nanmin(z_data))
        metrics.append(np.nanmean(disp_data))
        metrics.append(np.nanstd(disp_data))
        metrics.append(np.nanmedian(disp_data))
        metrics.append(np.percentile(disp_data, 80))
        metrics.append(np.percentile(disp_data, 90))
        metrics.append(np.percentile(disp_data, 95))
        metrics.append(np.nanmax(disp_data))
        metrics.append(np.nanmean(rigid_body_data))
        metrics.append(np.nanstd(rigid_body_data))
        metrics.append(np.nanmedian(rigid_body_data))
        metrics.append(np.nanmax(rigid_body_data))
        return metrics

    for fraction in range(len(fraction_xml_paths)):
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset = (
            read_motion_data(fraction_xml_paths[fraction])
        )
        time, potential_diff, rigid_body, x_offset, y_offset, z_offset, pauses = (
            modify_motion_data(
                time,
                potential_diff,
                rigid_body,
                x_offset,
                y_offset,
                z_offset,
                remove_duplicates=True,
                remove_variable_pauses=False,
            )
        )
        (
            zero_time,
            zero_potential_diff,
            zero_rigid_body,
            zero_x_offset,
            zero_y_offset,
            zero_z_offset,
            zero_pauses,
        ) = modify_motion_data(
            time,
            potential_diff,
            rigid_body,
            x_offset,
            y_offset,
            z_offset,
            remove_variable_pauses=False,
            zero_reference_point=True,
        )
        disp_time, disp = calculate_vector_displacements(
            zero_time, zero_x_offset, zero_y_offset, zero_z_offset
        )
        time_combined = [*time_combined, *zero_time]
        potential_diff_combined = [*potential_diff_combined, *potential_diff]
        rigid_body_combined = [*rigid_body_combined, *rigid_body]
        disp_time_combined = [*disp_time_combined, *disp_time]
        x_displacements = [*x_displacements, *zero_x_offset]
        y_displacements = [*y_displacements, *zero_y_offset]
        z_displacements = [*z_displacements, *zero_z_offset]
        v_displacements = [*v_displacements, *disp]
        fraction_results = [fraction + 1]
        fraction_results.append(time[-1] - time[0])
        fraction_results.append(len(pauses) + 1)
        pause_duration = 0
        for pause in pauses:
            pause_duration += pause[1] - pause[0]
        fraction_results.append(pause_duration)
        fraction_results.append((time[-1] - time[0]) - pause_duration)
        fraction_results = fraction_results + extract_metrics(
            zero_x_offset, zero_y_offset, zero_z_offset, disp, rigid_body
        )
        results.append(fraction_results)
    return_result = pd.DataFrame(results, columns=header)
    cumulative_result = [all_fraction_label]
    cumulative_result.append(return_result[header[1]].sum())
    cumulative_result.append(return_result[header[2]].sum())
    cumulative_result.append(return_result[header[3]].sum())
    cumulative_result.append(return_result[header[4]].sum())
    cumulative_result = cumulative_result + extract_metrics(
        x_displacements,
        y_displacements,
        z_displacements,
        v_displacements,
        rigid_body_combined,
    )
    return_result.loc[len(return_result.index)] = cumulative_result
    return return_result


def analyse_patient_motion_data(patient_path: str):
    """Analyse motion data contained within specific patient directory.

    Parameters
    ----------
    patient_path : str
        The path containing the *motionData.xml files to be plotted.

    Notes
    -----
    The specified directory should correspond with those cached in C:/tomo/da/pts/URnumber/
    when patient data is loaded within the Delivery Analysis tool.
    """
    motion_paths = patient_fraction_xml_lists(patient_path)
    return analyse_motion_data(motion_paths)
