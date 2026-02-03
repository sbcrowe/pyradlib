# -*- coding: utf-8 -*-
"""Motion analysis module.

This module provides functionality for processing of motion data from Synchrony treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import binascii
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import os
import pandas as pd
import pydicom
import xml.etree.ElementTree as et
from datetime import datetime
from scipy import ndimage

# define global variables
_mlc_thickness_at_iso = 6.25
_channel_thickness_at_iso = 1.24 * 850 / 1371.5


class Sinogram:

    def __init__(
        self,
        data: npt.ArrayLike,
        projections: npt.ArrayLike = None,
        fragment_start_projections: npt.ArrayLike = None,
        gantry_angles: npt.ArrayLike = None,
        times: npt.ArrayLike = None,
        fragment_start_datetimes: datetime = None,
        fragment_durations: npt.ArrayLike = None,
    ):
        """Initialises an object corresponding to a sinogram.

        Parameters
        ----------
        data : array_like
            Leaf open times or detector channel signal data.
        projections : array_like, optional
            Projection numbers for the sinogram data.
        fragment_start_projections : array_like, optional
            Projection numbers of the start of each beam delivery fragment.
        gantry_angles : array_like, optional
            Gantry angles corresponding to start of each projection.
        times : array_like, optional
            Times, in sec, corresponding to start of each projection,
            relative to start of initial fragment of session.
        fragment_start_datetimes : array_like, optional
            Date and time of the start of each beam delivery fragment.
        fragment_durations : array_like, optional
            The duration, in sec, of each beam delivery fragment.
        """
        self.data = data
        self.projections = projections
        self.fragment_start_projections = fragment_start_projections
        self.gantry_angles = gantry_angles
        self.times = times
        self.fragment_start_datetimes = fragment_start_datetimes
        self.fragment_durations = fragment_durations
        if self.projections is None:
            self.projections = np.array(range(len(data)))

    @classmethod
    def from_dcm(cls, path):
        """Reads planned sinogram from a DICOM RTPLAN file.

        Parameters
        ----------
        path : str
            Path to the .dcm file.

        Returns
        -------
        Sinogram
            The DICOM RTPLAN file is returned as a two-dimensional numpy array
            encapsulated with helper functions.
        """
        ds = pydicom.dcmread(path)
        gantry_angles = []
        sinogram_data = []
        for cp in ds.BeamSequence[0].ControlPointSequence:
            if (0x300D, 0x10A7) in cp:
                gantry_angles += [float(cp.GantryAngle)]
                sinogram_data += [
                    [
                        float(x)
                        for x in cp[0x300D, 0x10A7].value.decode("utf-8").split("\\")
                    ]
                ]
        return cls(np.array(sinogram_data), gantry_angles=np.array(gantry_angles))

    @classmethod
    def from_det(cls, path, clip: bool = True, flip: bool = True):
        """Reads detector sinogram from a .det file.

        Parameters
        ----------
        path : str
            Path to the .det file.
        clip : bool, optional
            Flag to indicate whether non-detector channels should be removed.
        flip : bool, optional
            Flag to indicate whether channels should be mirrored left-to-right,
            for ease of visual comparison against leaf sinograms.

        Returns
        -------
        Sinogram
            The det file is returned as a two-dimensional numpy array encapsulated
            with helper functions.

        Notes
        -----
        The function converts the binary-to-ASCII (hexadecimal) encoded data to a
        sinogram. This depends on the presence of fragment counts to determine how
        many projections exist for each leaf channel. For each leaf channel, the
        number of projections with open leaves is calculated from paired '>f8' tau
        values.
        """
        rows = int((os.path.getsize(path) - 16) / (640 * 4))
        data = np.fromfile(path, dtype=np.float32, offset=16).reshape((rows, 640))
        if clip:
            data = data[:, 0:576]
        if flip:
            data = np.flip(data, axis=1)
            if not clip:
                data = data.roll(data, -64, axis=1)
        return cls(data)

    @classmethod
    def from_dplan(cls, path: str, timing_path: str = None):
        """Reads planned or telemetry sinogram from a .dplan file.

        Parameters
        ----------
        path : str
            Path to the .dplan file.
        timing_path : str, optional
            Path to the telemetry timing .dat file.

        Returns
        -------
        Sinogram
            The dplan file is returned as a two-dimensional numpy array encapsulated
            with helper functions, with timing data if dat file is specified.

        Notes
        -----
        The function converts the binary-to-ASCII (hexadecimal) encoded data to a
        sinogram. This depends on the presence of fragment counts to determine how
        many projections exist for each leaf channel. For each leaf channel, the
        number of projections with open leaves is calculated from paired '>f8' tau
        values.
        """
        with open(path) as text_file:
            dplan_text = text_file.read()
        dplan_text = dplan_text.split("\n")
        seconds_per_tau = float(dplan_text[4].replace("fragment.secondsPerTau=", ""))
        start_tau = int(dplan_text[5].replace("fragment.startTau=", ""))
        end_tau = int(dplan_text[6].replace("fragment.endTau=", ""))
        tau_diff = end_tau - start_tau
        counts = [
            int(i) for i in dplan_text[7].replace("fragment.counts=", "").split(" ")
        ]
        dt = np.dtype(">f8")
        bin = np.frombuffer(binascii.unhexlify("".join(dplan_text[9:-12])), dtype=dt)
        data = np.zeros((tau_diff, 64))
        offset = 0
        for index, count in enumerate(counts[:-7]):
            if count > 0:
                # read bytes
                leaf_open_tau = bin[offset : offset + count]
                leaf_open_tau = leaf_open_tau - start_tau
                leaf_open_tau = leaf_open_tau.reshape(int(len(leaf_open_tau) / 2), 2)
                for tau in leaf_open_tau:
                    data[int(tau[0]) - 1, index] = tau[1] - tau[0]
                offset = offset + count
        projections = np.arange(start_tau, end_tau, 1)
        # calculate gantry, assuming projection 1 is centred on gantry 0
        gantry_angles = ((projections % 51) * (360 / 51) - (180 / 51)) % 360
        times = np.array(range(len(projections))) * seconds_per_tau
        if timing_path is not None:
            with open(timing_path) as timing_path_text_file:
                timing_text = timing_path_text_file.read()
            timing_text = timing_text.split("\n")
            fragment_start_timestamps = np.array(
                [
                    x + y / 1000000
                    for x, y in zip(
                        list(map(int, timing_text[3].split(","))),
                        list(map(int, timing_text[4].split(","))),
                    )
                    if x > 0
                ]
            )
            fragment_start_datetimes = [
                datetime.fromtimestamp(x) for x in fragment_start_timestamps
            ]
            fragment_breaks = np.array(
                [x for x in list(map(float, timing_text[5].split(","))) if x > 0]
            )
            fragment_durations = np.array(
                [
                    (y - x) * seconds_per_tau
                    for x, y in zip(
                        list(map(float, timing_text[5].split(","))),
                        list(map(float, timing_text[6].split(","))),
                    )
                    if x > 0
                ]
            )
            times += np.piecewise(
                projections,
                [projections >= fragment_break for fragment_break in fragment_breaks],
                fragment_start_timestamps
                - fragment_start_timestamps[0]
                - (fragment_breaks - fragment_breaks[0]) * seconds_per_tau,
            )
        return cls(
            data,
            projections=projections,
            fragment_start_projections=fragment_breaks,
            gantry_angles=gantry_angles,
            times=times,
            fragment_start_datetimes=fragment_start_datetimes,
            fragment_durations=fragment_durations,
        )

    @classmethod
    def from_bin(cls, path: str, channels: int = 64):
        """Reads planned, telemetry, difference or detector sinogram from a .bin
        file exported from the Delivery Analysis software.

        Parameters
        ----------
        path : str
            Path to the .bin file.

        Returns
        -------
        Sinogram
            The csv file is returned as a two-dimensional numpy array encapsulated
            with helper functions.

        Notes
        -----
        The format of the bin file is an int32 describing number of projections, an
        int32 describing number of channels (e.g., 64 for leaf sinogram or 640 for
        detector sinogram), an int32 with unknown use ("1"), then an array of
        float32 values with dimensions (num_projections, 4) including gantry angles
        and couch positions, then an array of float32 values with dimensions (
        num_projections, num_channels) including leaf open times or detector counts.
        """
        num_projections, num_channels = np.fromfile(path, dtype=np.int32, count=2)
        header = np.fromfile(
            path, dtype=np.float32, offset=12, count=num_projections * 4
        ).reshape((num_projections, 4))
        gantry_angles = (header[:, 0] - (180 / 51)) % 360
        couch_positions = header[:, 3]
        data = np.fromfile(
            path, dtype=np.float32, offset=12 + num_projections * 4 * 4
        ).reshape((num_projections, num_channels))
        if data.shape[1] == 640:
            data = data[:, 0:576]
            data = np.flip(data, axis=1)
        return cls(data, gantry_angles=gantry_angles)

    @classmethod
    def from_csv(cls, path: str):
        """Reads planned, telemetry, difference or detector sinogram from a .csv
        file exported from the Delivery Analysis software.

        Parameters
        ----------
        path : str
            Path to the .csv file.

        Returns
        -------
        Sinogram
            The csv file is returned as a two-dimensional numpy array encapsulated
            with helper functions.
        """
        data = np.loadtxt(path, delimiter=",", skiprows=1)[:, 3:]
        if data.shape[1] == 640:
            data = data[:, 0:576]
            data = np.flip(data, axis=1)
        projections = np.loadtxt(path, delimiter=",", skiprows=1)[:, 0]
        gantry_angles = (
            np.loadtxt(path, delimiter=",", skiprows=1)[:, 1] - (180 / 51)
        ) % 360
        return cls(data, projections=projections, gantry_angles=gantry_angles)

    def __add__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError("Unsupported operand type for +")
        return type(self)(self.data + other.data, projections=np.copy(self.projections))

    def __len__(self):
        return len(self.data)

    def __sub__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError("Unsupported operand type for -")
        return type(self)(self.data - other.data, projections=np.copy(self.projections))

    def adapt_to_motion(self, motion):
        # TODO implement this function
        raise NotImplementedError("Motion adaptation in development")
        adapted_data = np.copy(self.data)
        # define gantry angles matching motion
        # calculate beams eye view motion
        # for all points in sinogram, calculate shift
        return Sinogram(
            adapted_data,
            projections=np.copy(self.projections),
            gantry_angles=np.copy(self.gantry_angles),
        )

    def metrics(self):
        header_row = ["Projections", "totalLOT", "TA"]
        treatment_area_leaves = 0
        for projection in self.data:
            if np.sum(projection) > 0:
                left_index_open = next((i for i, x in enumerate(projection) if x), None)
                right_index_open = 63 - next(
                    (i for i, x in enumerate(reversed(projection)) if x), None
                )
                treatment_area_leaves += right_index_open - left_index_open + 1
        return [
            self.data.shape[0],
            np.sum(self.data) / self.data.shape[0],
            treatment_area_leaves / self.data.shape[0],
        ]

    def plot_sinogram(self, ax, title: str = None, vmin: float = 0, vmax: float = 1):
        """Plots sinogram data to Matplotlib specified axis."""
        out = ax.imshow(
            self.data,
            vmin=0,
            vmax=1,
            extent=[-200, 200, self.projections[-1], self.projections[0]],
        )
        if not title is None:
            ax.set_title(title)
        return out

    def plot_projection_lot(self, ax, title: str = None):
        out = ax.plot(np.sum(self.data, axis=1))
        return out

    def scale_to_isocentre(self):
        """Scales sinogram to 1 mm resolution at isocentre.

        Returns
        -------
        Sinogram
            The sinogram rescaled to 1 mm resolution at isocentre.

        Notes
        -----
        The function is intended to allow pixel-by-pixel comparison of planned or
        telemetry leaf open time and detector channel signals. This necessarily
        requires resampling of the data (e.g., each leaf will become 6.25 array
        elements). Due to the resampling, the values of each array element will
        not accurately reflect delivery parameters.
        """
        if len(self.data[1]) == 64:
            return Sinogram(
                ndimage.zoom(self.data, (1, _mlc_thickness_at_iso)), self.projections
            )
        elif len(self.data[1]) == 576:
            scaled_data = ndimage.zoom(self.data, (1, _channel_thickness_at_iso))
            shift_offset = (len(scaled_data[0]) - 400) / 2
            scaled_data = ndimage.shift(scaled_data, (0, -shift_offset))[:, 0:400]
            return Sinogram(scaled_data, self.projections)
        return None

    def to_csv(self, path):
        """Writes sinogram leaf open times or channel signal data to CSV file.

        Parameters
        ----------
        path : str
            Path of the CSV file to be written.
        """
        if self.data.shape[1] <= 64:
            header_row = ["Projection"] + [
                "Leaf %d LOT" % (i + 1) for i in range(self.data.shape[1])
            ]
        else:
            header_row = ["Projection"] + [
                "Channel %d count" % (i + 1) for i in range(self.data.shape[1])
            ]
        projection_column = np.arange(1, self.data.shape[0] + 1).reshape((-1, 1))
        np.savetxt(
            path,
            np.concatenate([projection_column, self.data], axis=1),
            fmt="%i" + (",%8f" * self.data.shape[1]),
            delimiter=",",
            header=",".join(header_row),
        )

    def to_dcm(self, path, template_path):
        """Writes sinogram leaf open times to DICOM RTPLAN file, using an existing
        DICOM RTPLAN file as a template.

        Parameters
        ----------
        path : str
            Path of the DICOM RTPLAN file to be written.
        template_path : str
            Path of the DICOM RTPLAN file to be used as a template.

        Notes
        -----
        Depending on how the DICOM RTPLAN has been collected, e.g., either exported
        from the TPS or collected from a PDE or delivery analysis cache, the first
        control point may or may not have associated leaf open times.
        """
        ds = pydicom.dcmread(template_path)
        if len(self.data) == len(ds.BeamSequence[0].ControlPointSequence) - 1:
            for cp, lot in zip(
                ds.BeamSequence[0].ControlPointSequence[0:-1], self.data
            ):
                cp[0x300D, 0x10A7].value = "//".join([str(x) for x in lot]).encode(
                    "utf-8"
                )
        elif len(self.data) == len(ds.BeamSequence[0].ControlPointSequence) - 2:
            for cp, lot in zip(
                ds.BeamSequence[0].ControlPointSequence[1:-1], self.data
            ):
                cp[0x300D, 0x10A7].value = "//".join([str(x) for x in lot]).encode(
                    "utf-8"
                )
        else:
            raise RuntimeError(
                "Template DICOM RTPLAN has incorrect number of control points"
            )
        ds.save_as(path)

    def to_png(self, path, title=None):
        """Writes sinogram leaf open times or channel signal data to PNG file.

        Parameters
        ----------
        path : str
            Path of the PNG file to be written.
        """
        fig, ax = plt.subplots(1, 1, figsize=(4.5, 3 * self.data.shape[0] / 400))
        im = ax.imshow(
            self.data, extent=[-200, 200, self.projections[-1], self.projections[0]]
        )
        cbar = plt.colorbar(im, ax=ax)
        if self.data.shape[1] <= 64:
            cbar_title = "LOT"
        else:
            cbar_title = "Counts"
        cbar.ax.set_title(cbar_title)
        if title is not None:
            fig.suptitle(title)
        fig.supxlabel("Distance from isocentre (mm)")
        fig.supylabel("Projection")
        fig.tight_layout(pad=1)
        fig.savefig(path)
        plt.close(fig)


def plot_sinograms(
    path: str,
    plan_sinogram: Sinogram,
    telemetry_sinogram: Sinogram = None,
    detector_sinogram: Sinogram = None,
    difference_scale: list = [1],
    title: str = None,
):
    num_plots = 1
    if telemetry_sinogram is not None:
        num_plots = num_plots + 1
        if (plan_sinogram.data.shape == telemetry_sinogram.data.shape) and len(
            difference_scale
        ) > 0:
            num_plots = num_plots + len(difference_scale)
    if detector_sinogram is not None:
        num_plots = num_plots + 1
    fig_width = 4.5 * num_plots
    fig_height = 3 * plan_sinogram.data.shape[0] / 400
    fig, axes = plt.subplots(1, num_plots, figsize=(fig_width, fig_height))
    im_plan = axes[0].imshow(
        plan_sinogram.data,
        vmin=0,
        vmax=1,
        extent=[-200, 200, plan_sinogram.data.shape[0], 0],
    )
    axes[0].set_title("Plan")
    cbar_plan = plt.colorbar(im_plan, ax=axes[0])
    cbar_plan.ax.set_title("LOT")
    next_fig = 1
    if telemetry_sinogram is not None:
        im_telem = axes[1].imshow(
            telemetry_sinogram.data,
            vmin=0,
            vmax=1,
            extent=[-200, 200, telemetry_sinogram.data.shape[0], 0],
        )
        axes[next_fig].set_title("Telemetry")
        cbar_telem = plt.colorbar(im_telem, ax=axes[next_fig])
        cbar_telem.ax.set_title("LOT")
        next_fig = next_fig + 1
    if (plan_sinogram.data.shape == telemetry_sinogram.data.shape) and len(
        difference_scale
    ) > 0:
        for scale in difference_scale:
            im_diff = axes[next_fig].imshow(
                telemetry_sinogram.data - plan_sinogram.data,
                vmin=-1 / scale,
                vmax=1 / scale,
                extent=[-200, 200, plan_sinogram.data.shape[0], 0],
            )
            axes[next_fig].set_title("Difference")
            cbar_diff = plt.colorbar(im_diff, ax=axes[next_fig])
            cbar_diff.ax.set_title("ΔLOT")
            next_fig = next_fig + 1
    if detector_sinogram is not None:
        im_detec = axes[next_fig].imshow(
            detector_sinogram.data,
            vmin=0,
            extent=[-221.33, 221.33, detector_sinogram.data.shape[0], 0],
        )
        axes[next_fig].set_title("Detector")
        cbar_detec = plt.colorbar(im_detec, ax=axes[next_fig])
        cbar_detec.ax.set_title("Counts")
    if title is not None:
        fig.suptitle(title)
    fig.supxlabel("Distance from isocentre (mm)")
    fig.supylabel("Projection")
    fig.tight_layout(pad=1)
    fig.savefig(path)
    plt.close(fig)


def plot_projection_leaf_open_time(
    path: str,
    plan_sinogram: Sinogram,
    telemetry_sinograms: npt.ArrayLike = [],
    title: str = None,
    number_columns: int = 2,
    wide_plot: bool = True,
):
    if len(telemetry_sinograms) > 0:
        # plot sinograms
        ncols = np.min([number_columns, len(telemetry_sinograms)])
        nrows = int((len(telemetry_sinograms) + (ncols - 1)) / ncols)
        if wide_plot:
            width_multiplier = 8
        else:
            width_multiplier = 2
        fig, axs = plt.subplots(
            ncols=ncols,
            nrows=nrows,
            sharex=True,
            sharey=True,
            gridspec_kw={"hspace": 0, "wspace": 0},
            constrained_layout=True,
            figsize=(ncols * width_multiplier + 1, nrows * 2 + 1),
        )
        for i in range(len(telemetry_sinograms)):
            telemetry_sinogram = telemetry_sinograms[i]
            col = int((i) % ncols)
            row = int((i) / ncols)
            if nrows > 1:
                (l1,) = axs[row, col].plot(np.sum(plan_sinogram.data, axis=1))
                (l2,) = axs[row, col].plot(np.sum(telemetry_sinogram.data, axis=1))
            elif ncols > 1:
                (l1,) = axs[col].plot(np.sum(plan_sinogram.data, axis=1))
                (l2,) = axs[col].plot(np.sum(telemetry_sinogram.data, axis=1))
            else:
                (l1,) = axs.plot(np.sum(plan_sinogram.data, axis=1))
                (l2,) = axs.plot(np.sum(telemetry_sinogram.data, axis=1))
            plotted = [l1, l2]
            legend = ["Plan", "Telemetry"]
        fig.legend(
            plotted,
            legend,
            loc="center right",
            ncol=ncols,
            fancybox=True,
            bbox_to_anchor=(1.05, 1),
        )
        fig.supxlabel("Projection")
        fig.supylabel("Integrated leaf open time (s)")
        if title is not None:
            fig.suptitle(title)
        fig.savefig(path, bbox_inches="tight")
        plt.close(fig)
