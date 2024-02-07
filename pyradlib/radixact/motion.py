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
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import os
import re
import xml.etree.ElementTree as et


def read_motion_data(xml_filepaths: npt.ArrayLike):
    """Extracts motion data from *motionData.xml files produced for the Radixact Delivery Analysis tool.

    Parameters
    ----------
    xml_filepaths : array_like
        List of paths to XML files containing motion data.

    Returns
    -------
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset : array_like
        Output arrays, containing paired timestamps, potential differences, rigid body differences, and IEC-X, IEC-Y and IEC-Z offsets (in sec, mm).

    Notes
    -----
    These files are cached in C:/tomo/da/pts/URnumber/*motionData.xml when patient data is loaded within the Delivery Analysis tool. There may be multiple files per fraction.

    Pauses in fractions (e.g. due to user, Synchrony, or system errors) are indicated by numpy.nan. This allows more accurate plotting of the data using matplotlib.
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
    """Converts motion data from *motionData.xml files produced for the Radixact Delivery Analysis tool to CSV files.

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
):
    """Modify extracted motion data for improved presentation.

    Parameters
    ----------
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset : array_like
        Input arrays, containing paired timestamps, potential differences, rigid body differences, and IEC-X, IEC-Y and IEC-Z offsets (in sec, mm).
    remove_duplicates : bool, optional
        Flag indicating removal of duplicate data points, including at start of treatment, for reduced data density and profile length (default is True).
    remove_variable_pauses : bool, optional
        Flag indicating replacement of variable length pauses in treatment with short fixed length, defined by pause_offset_time (default is True).
    replacement_pause_length : float, optional
        Time period used when replacing variable pauses in beam delivery, e.g., 5 sec (default).

    Returns
    -------
    time, potential_diff, rigid_body, x_offset, y_offset, z_offset : array_like
        Output arrays, containing adjusted paired timestamps, potential differences, rigid body differences, and IEC-X, IEC-Y and IEC-Z offsets (in sec, mm).
    pause_times : array_like
        Output array, containing 0 or more paired timestamps indicating when pauses in treatment occur, for plotting (in sec, mm).
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
    return (
        adjusted_time,
        new_potential_diff,
        new_rigid_body,
        new_x_offset,
        new_y_offset,
        new_z_offset,
        pauses,
    )


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
):
    """ Plot collection of motion data spanning multiple fractions.
    
    Parameters
    ----------
    fraction_xml_paths : array_like
        List containing lists of *motionData.xml paths for each fraction of the treatment.
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
        Flag indicating removal of duplicate data points, including at start of treatment, for reduced data density and profile length (default is True).
    remove_variable_pauses : bool, optional
        Flag indicating replacement of variable length pauses in treatment with short fixed length, defined by pause_offset_time (default is True).
    replacement_pause_length : float, optional
        Time period used when replacing variable pauses in beam delivery, e.g., 5 sec (default).
    pause_colour : str, optional
        Color to use to indicate pauses in treatment, e.g., 'lightgrey' (default).
    """
    nrows = int((len(fraction_xml_paths) + 4) / number_columns)
    fig, axs = plt.subplots(
        ncols=number_columns,
        nrows=nrows,
        sharex=share_x_axis,
        sharey=True,
        gridspec_kw={"hspace": 0, "wspace": 0},
        constrained_layout=True,
        figsize=(11, 8),
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
            )
        )
        col = int(fraction % number_columns)
        row = int(fraction / number_columns)
        (l1,) = axs[row, col].plot(time, x_offset)
        (l2,) = axs[row, col].plot(time, y_offset)
        (l3,) = axs[row, col].plot(time, z_offset)
        plotted = [l1, l2, l3]
        legend = ["IEC-X", "IEC-Y", "IEC-Z"]
        if plot_potential_diff:
            (l4,) = axs[row, col].plot(time, potential_diff)
            plotted.append(l4)
            legend.append("Potential Diff")
        if plot_rigid_body:
            (l5,) = axs[row, col].plot(time, rigid_body)
            plotted.append(l5)
            legend.append("Rigid Body")
        for pause_time in pauses:
            axs[row, col].axvspan(pause_time[0], pause_time[1], color=pause_color, lw=0)
        axs[row, col].set_title(fraction + 1, y=1.0, pad=-14)
    # plt.ylim(-_plot_maximum_motion,_plot_maximum_motion)
    fig.legend(
        plotted,
        legend,
        loc="lower center",
        ncol=5,
        fancybox=True,
        shadow=True,
        bbox_to_anchor=(0, -0.05, 1, 1),
    )
    fig.supxlabel("Time (s)")
    fig.supylabel("Distance (mm)")
    # remove empty grids
    if len(fraction_xml_paths) % number_columns > 0:
        axs[-1, -1].set_axis_off()
        axs[-2, -1].xaxis.set_tick_params(which="both", labelbottom=True)
        if len(fraction_xml_paths) % number_columns < 4:
            axs[-1, -2].set_axis_off()
            axs[-2, -2].xaxis.set_tick_params(which="both", labelbottom=True)
            if len(fraction_xml_paths) % number_columns < 3:
                axs[-1, -3].set_axis_off()
                axs[-2, -3].xaxis.set_tick_params(which="both", labelbottom=True)
                if len(fraction_xml_paths) % number_columns < 2:
                    axs[-1, -4].set_axis_off()
                    axs[-2, -4].xaxis.set_tick_params(which="both", labelbottom=True)
    fig.suptitle(title)
    plt.savefig(png_filepath)


def plot_patient_data(patient_path: str, png_path: str, title: str):
    """ Plot motion data contained within specific patient directory.

    Parameters
    ----------
    patient_path : str
        The path containing the *motionData.xml files to be plotted.
    png_filepath : str
        Path for PNG file to be written.
    title : str
        Title text to use for figure.

    Notes
    -----
    The specified directory should correspond with those cached in C:/tomo/da/pts/URnumber/ when patient data is loaded within the Delivery Analysis tool.
    """
    motion_files = sorted(
        glob.glob(os.path.join(patient_path, "*motionData.xml")),
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
    plot_motion_data(list(motion_paths.values()), png_path, title)
