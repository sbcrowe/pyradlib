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
import os
from functools import cached_property

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import polars as pl
import pydicom
from matplotlib import ticker
from scipy import stats
from scipy.ndimage import center_of_mass, shift

from pyradlib.radixact.motion import RadixactSynchronyMotion
from pyradlib.radixact.timing import RadixactTiming


class RadixactSinogram:
    # region Constructors

    def __init__(self, data: pl.DataFrame, uid: str = None):
        """Initialises an object corresponding to a sinogram.

        Parameters
        ----------
        data : DataFrame
            DataFrame containing projections, and optionally, gantry angles,
            datetimes, leaf open times, and detector channel signal data.
        uid : str, optional
            Unique identifier associated with the data. Default is None.
        """
        self._df = data
        self._uid = uid

    @classmethod
    def from_bin(cls, path: str | os.PathLike, uid: str = None) -> RadixactSinogram:
        """Reads planned, telemetry, difference or detector sinogram from a .bin
        file exported from the Delivery Analysis software.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the .bin file.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

        Returns
        -------
        RadixactSinogram
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
        df = pl.DataFrame(
            {
                "projection": np.arange(num_projections) + 1,
                "gantry_angle": gantry_angles,
                "couch_position": couch_positions,
                "channel_values": data,
            }
        )
        return cls(df, uid)

    @classmethod
    def from_dcm(cls, path: str | os.PathLike, uid: str = None) -> RadixactSinogram:
        """Reads planned sinogram from a DICOM RTPLAN file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the .dcm file.
        uid : str, optional
            Unique identifier associated with the data. Default is None, in which case
            a UID from the DICOM file will be used.

        Returns
        -------
        RadixactSinogram
            The DICOM RTPLAN file is returned as a two-dimensional numpy array
            encapsulated with helper functions.
        """
        ds = pydicom.dcmread(path)
        data = []
        for cp in ds.BeamSequence[0].ControlPointSequence:
            if (0x300D, 0x10A7) in cp:
                data.append(
                    {
                        "gantry_angle": float(cp.GantryAngle),
                        "channel_values": [
                            float(x)
                            for x in cp[0x300D, 0x10A7]
                            .value.decode("utf-8")
                            .split("\\")
                        ],
                    }
                )
        df = pl.DataFrame(data)
        if uid is None:
            # TODO Extract a UID from DCM file.
            dicom_uid = None
            return cls(df, dicom_uid)
        else:
            return cls(df, uid)

    @classmethod
    def from_delivery_analysis_csv(
        cls, path: str | os.PathLike, uid: str = None
    ) -> RadixactSinogram:
        """Reads planned, telemetry, difference or detector sinogram from a .csv
        file exported from the Delivery Analysis software.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the .csv file.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

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
        df = pl.DataFrame(
            {
                "projection": projections,
                "gantry_angle": gantry_angles,
                "channel_values": data,
            }
        )
        return cls(df)

    @classmethod
    def from_det(cls, path: str | os.PathLike, uid: str = None) -> RadixactSinogram:
        """Reads detector sinogram from a .det file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the .det file.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

        Returns
        -------
        RadixactSinogram
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
        data = data[:, 0:576]
        data = np.flip(data, axis=1)
        df = pl.DataFrame(
            {
                "projection": np.arange(rows) + 1,
                "channel_values": data,
            }
        )
        return cls(df, uid)

    @classmethod
    def from_dplan(cls, path: str | os.PathLike, uid: str = None) -> RadixactSinogram:
        """Reads planned or telemetry sinogram from a .dplan file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the .dplan file.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

        Returns
        -------
        RadixactSinogram
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
        # Read dplan file
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
        # Calculate leaf open times
        for index, count in enumerate(counts[:-7]):
            if count > 0:
                leaf_open_tau = bin[offset : offset + count]
                leaf_open_tau = leaf_open_tau - start_tau
                leaf_open_tau = leaf_open_tau.reshape(int(len(leaf_open_tau) / 2), 2)
                for tau in leaf_open_tau:
                    data[int(tau[0]) - 1, index] = tau[1] - tau[0]
                offset = offset + count
        tau = np.arange(start_tau, end_tau, 1)
        # Calculate gantry angle, assuming projection 1 is centred on gantry 0
        gantry_angles = ((tau % 51) * (360 / 51) - (180 / 51)) % 360
        # Calculate nominal delivery times
        times = np.array(range(len(tau))) * seconds_per_tau
        # Return data if no telemetry timing available
        df = pl.DataFrame(
            {
                "projection": tau - tau[0] + 1,
                "tau": tau,
                "time": times,
                "gantry_angle": gantry_angles,
                "channel_values": data,
            }
        )
        return cls(df, uid)

    @classmethod
    def from_parquet(cls, path: str | os.PathLike, uid: str = None) -> RadixactSinogram:
        """Reads sinogram from a parquet file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the parquet file.
        uid : str, optional
            Unique identifier associated with the data. Default is None.

        Returns
        -------
        RadixactSinogram
            Two-dimensional numpy array encapsulated with helper functions.
        """
        df = pl.read_parquet(path)
        return cls(df, uid)

    # endregion

    # region Magic methods

    def __len__(self) -> int:
        """Returns number of sinogram projections.

        Returns
        -------
        int
            Number of projections.
        """
        return len(self._df)

    # endregion

    # region Properties

    @cached_property
    def metrics(self) -> pl.DataFrame:
        """Returns sinogram metrics.

        Returns
        -------
        pl.DataFrame
            The dataframe containing sinogram metrics.

        Notes
        -------
        The implementation was informed by Cavinato et al.'s "Quantitative assessment of
        helical tomotherapy plans complexity", DOI:10.1002//acm2.13781
        """
        # TODO Handle upsampled or detector sinograms.
        metrics_dict = {}
        projection_time = (
            self._df.get_column("time").to_numpy()[1]
            - self._df.get_column("time").to_numpy()[0]
        ) * 1000
        leaf_open_fractions = self.fractional_leaf_open_times()
        leaf_open_times = self.leaf_open_times()
        non_zero_leaf_open_times = leaf_open_times[leaf_open_times > 0]
        non_zero_leaf_open_fractions = leaf_open_fractions[leaf_open_fractions > 0]
        # Calculate projection statistics
        metrics_dict["projection_time"] = projection_time
        metrics_dict["num_projections"] = len(leaf_open_fractions)
        metrics_dict["num_rotations"] = len(leaf_open_fractions) / 51
        # Calculate absolute LOT statistics
        metrics_dict["num_non_zero_lots"] = len(non_zero_leaf_open_times)
        metrics_dict["mean_lot"] = np.mean(non_zero_leaf_open_times)
        metrics_dict["median_lot"] = np.median(non_zero_leaf_open_times)
        metrics_dict["mode_lot"] = stats.mode(non_zero_leaf_open_times)[0]
        metrics_dict["standard_deviation_lot"] = np.std(non_zero_leaf_open_times)
        metrics_dict["kurtosis_lot"] = stats.kurtosis(non_zero_leaf_open_times)
        metrics_dict["skewness_lot"] = stats.skew(non_zero_leaf_open_times)
        metrics_dict["min_lot"] = np.min(non_zero_leaf_open_times)
        metrics_dict["max_lot"] = np.max(non_zero_leaf_open_times)
        # Calculate cumulative absolute LOT number scores
        for clns_time in [100, 50, 30, 20]:
            metrics_dict["clns_" + str(clns_time)] = len(
                non_zero_leaf_open_times[non_zero_leaf_open_times < clns_time]
            ) / len(non_zero_leaf_open_times)
        for clns_pt_time in [50, 20]:
            metrics_dict["clns_pt_" + str(clns_pt_time)] = len(
                non_zero_leaf_open_times[
                    non_zero_leaf_open_times > (projection_time - clns_pt_time)
                ]
            ) / len(non_zero_leaf_open_times)
        # Calculate relative LOT statistics
        metrics_dict["mean_flot"] = np.mean(non_zero_leaf_open_fractions)
        metrics_dict["median_flot"] = np.median(non_zero_leaf_open_fractions)
        metrics_dict["mode_flot"] = stats.mode(non_zero_leaf_open_fractions)[0]
        metrics_dict["standard_deviation_flot"] = np.std(non_zero_leaf_open_fractions)
        metrics_dict["min_flot"] = np.min(non_zero_leaf_open_fractions)
        metrics_dict["max_flot"] = np.max(non_zero_leaf_open_fractions)
        # Calculate cumulative fractional lot number scores
        for cfns_fraction in [5, 10, 50, 75, 90]:
            metrics_dict["cfns_" + str(cfns_fraction)] = len(
                non_zero_leaf_open_fractions[
                    non_zero_leaf_open_fractions < cfns_fraction
                ]
            ) / len(non_zero_leaf_open_fractions)
        # TODO Calculate sinogram geometry statistics
        # L0NS, L1NS, CLS, CLSin, L2NS, nCC, lengthCC, TA, fDISC, CLSin,area, CLSin,disc
        # CLS in,area,disc, centroid
        treatment_area_leaves = 0
        for projection in leaf_open_fractions:
            if np.sum(projection) > 0:
                left_index_open = next((i for i, x in enumerate(projection) if x), None)
                right_index_open = 63 - next(
                    (i for i, x in enumerate(reversed(projection)) if x), None
                )
                treatment_area_leaves += right_index_open - left_index_open + 1
        treatment_area = treatment_area_leaves / len(leaf_open_fractions)
        metrics_dict["treatment_area"] = treatment_area
        # TODO Calculation sinogram modulation statistics
        # PSTV, LOTV, MI, nOC, EPSTV1,1, EPSTV0,1, EPSTV1,0, ELOTV1, ELOTV2, ELOTV3,
        # ELOTV4, ELOTV5, mSI, mdSI, sdSI, MSA
        return pl.from_dict(metrics_dict)

    # endregion

    # region Public methods

    def apply_motion_adaptation(
        self,
        timing: RadixactTiming,
        motion: RadixactSynchronyMotion,
        leaf_adaptation_mode: str = "threshold",
        leaf_adaptation_threshold: float = 4.3 / 6.25,
    ) -> RadixactSinogram:
        """Adapts sinogram to specified motion.

        Parameters
        ----------
        timing : RadixactTiming
            The telemetry timing to use for motion correlation.
        motion : RadixactSynchronyMotion
            The motion data to use to calculate expected adaptation.
        leaf_adaptation_mode : {'threshold', 'optimal'}, optional
            The mode parameter determines how the adaptation is calculated.
            'threshold'
                Adaptation occurs when target offset relative to leaf thickness in the
                beam's eye view exceeds tolerance defined by tval.
            'optimal'
                Adaptation is proportional to the target offset in the beam's eye view.
                Non-integer leaf adjustments are modelled.
        leaf_adaptation_threshold : float, optional
            The change in target offset relative to leaf thickness in the beam's eye
            view sufficient to trigger leaf adaptation. Default is 0.688 of a leaf
            thickness, corresponding to 4.3 mm, reported by Miura et al (2025).

        Returns
        -------
        RadixactSinogram
            The expected sinogram after adaptation to motion.
        """
        # Add timing data to sinogram.
        df = self.apply_timing(timing)._df
        time_interval = df["timestamp"][1] - df["timestamp"][0]
        df = df.with_columns(
            ((pl.col("timestamp") + time_interval / 2) / 1000).alias(
                "adaptation_timestamp"
            )
        )
        # Interpolate motion
        df = df.with_columns(
            target_offset_x=pl.Series(
                np.interp(
                    x=df["adaptation_timestamp"],
                    xp=motion._df["timestamp"],
                    fp=motion._df["target_offset_x"],
                )
            )
        )
        df = df.with_columns(
            target_offset_y=pl.Series(
                np.interp(
                    x=df["adaptation_timestamp"],
                    xp=motion._df["timestamp"],
                    fp=motion._df["target_offset_y"],
                )
            )
        )
        df = df.with_columns(
            target_offset_z=pl.Series(
                np.interp(
                    x=df["adaptation_timestamp"],
                    xp=motion._df["timestamp"],
                    fp=motion._df["target_offset_z"],
                )
            )
        )
        # Calculate offset in beams eye view
        gantry_interval = (df["gantry_angle"][1] - df["gantry_angle"][0]) % 360
        # df = df.with_columns(
        #    adaptation_gantry_angle = (((df["projection"]) * gantry_interval)/2 + df["gantry_angle"][0]) % 360
        # )
        df = df.with_columns(
            ((pl.col("gantry_angle") + gantry_interval / 2) % 360).alias(
                "adaptation_gantry_angle"
            )
        )
        df = df.with_columns(
            target_offset_beams_eye_view=np.cos(
                np.deg2rad(df["adaptation_gantry_angle"])
            )
            * df["target_offset_x"]
            - np.sin(np.deg2rad(df["adaptation_gantry_angle"])) * df["target_offset_z"]
        )
        df = df.with_columns(float_leaf_shift=df["target_offset_beams_eye_view"] / 6.25)
        # Apply adaptation
        if leaf_adaptation_mode == "optimal":
            df = df.with_columns(pl.col("float_leaf_shift").alias("leaf_shift"))
            optimal_leaf_open_times = np.zeros_like(df["channel_values"])
            for i, s in enumerate(df["leaf_shift"]):
                optimal_leaf_open_times[i] = shift(
                    df["channel_values"][i], s, mode="constant", cval=0
                )
            df = df.with_columns(channel_values=optimal_leaf_open_times)
        elif leaf_adaptation_mode == "threshold":
            predicted_leaf_shift = []
            curr_shift = 0
            for row in df.iter_rows(named=True):
                if (
                    abs(row["float_leaf_shift"] - curr_shift)
                    >= leaf_adaptation_threshold
                ):
                    new_shift = int(np.rint(row["float_leaf_shift"]))
                    predicted_leaf_shift.append(new_shift)
                    curr_shift = new_shift
                else:
                    predicted_leaf_shift.append(curr_shift)
            df = df.with_columns(pl.Series("leaf_shift", predicted_leaf_shift))
            predicted_leaf_open_times = np.zeros_like(df["channel_values"])
            for i, s in enumerate(df["leaf_shift"]):
                if s > 0:
                    predicted_leaf_open_times[i][s:] = df["channel_values"][i][:-s]
                elif s < 0:
                    predicted_leaf_open_times[i][:s] = df["channel_values"][i][-s:]
                else:
                    predicted_leaf_open_times[i] = df["channel_values"][i]
            df = df.with_columns(channel_values=predicted_leaf_open_times)
        return type(self)(df)

    def apply_timing(self, timing: RadixactTiming) -> RadixactSinogram:
        """Associate telemetry timing information with sinogram data.

        Parameters
        ----------
        timing : RadixactTiming
            Telemetry timing data.

        Returns
        -------
        RadixactSinogram
            Sinogram with associated telemetry timing information.
        """
        df = self._df
        df = df.with_columns(
            timestamp=pl.Series(
                np.interp(x=df["tau"], xp=timing._df["tau"], fp=timing._df["timestamp"])
            )
        )
        df = df.with_columns(
            pl.from_epoch("timestamp", time_unit="us").alias("datetime")
        )
        return type(self)(df)

    def closed_leaf_projections(self, flot_open_threshold=0.1) -> npt.ArrayLike:
        """Returns sinogram row mask for projections with all leaves closed.

        Parameters
        ----------
        flot_open_threshold : float, optional
            Threshold to use to determine if any leaves are open across projection.
            Default 0.2

        Returns
        -------
        npt.ArrayLike
            Sinogram row mask containing True and False values, for all leaves being
            closed within a given projection.
        """
        result = []
        self_data = self._df["channel_values"].to_numpy()
        for self_row in self_data:
            if np.sum(self_row) < flot_open_threshold:
                result.append(True)
            else:
                result.append(False)
        return result

    def downsample(self) -> RadixactSinogram:
        """Downsamples sinogram with fractional projctions to whole projections.

        Returns
        -------
        RadixactSinogram
            The sinogram downsampled to whole projections.
        """
        # Calculate new projection, tau, time and gantry angles.
        new_projection = self._df["projection"].unique()
        new_tau = self._df["tau"].cast(pl.Int64).unique()
        new_time_interval = (self._df["time"][1] - self._df["time"][0]) / (
            self._df["tau"][1] - self._df["tau"][0]
        )
        new_time = (new_projection.to_numpy() - 1) * new_time_interval
        new_gantry_angles = ((new_tau % 51) * (360 / 51) - (180 / 51)) % 360
        # Create new data frame for updated leaf open times.
        df = pl.DataFrame(
            {
                "projection": new_projection,
                "tau": new_tau,
                "time": new_time,
                "gantry_angle": new_gantry_angles,
            }
        )
        # Calculate new leaf open times for downsampled projections.
        new_leaf_open_times = []
        for downsampled_row in df.iter_rows(named=True):
            curr_rows = self._df.filter(
                pl.col("projection") == downsampled_row["projection"]
            )
            mean_leaf_open_times = np.mean(
                curr_rows["channel_values"].to_numpy(), axis=0
            )
            new_leaf_open_times.append(mean_leaf_open_times)
        df = df.with_columns(channel_values=pl.Series(np.array(new_leaf_open_times)))
        return type(self)(df)

    def estimated_leaf_adaptation(
        self, other: RadixactSinogram, mode="first-open", flot_open_threshold=0.2
    ) -> npt.ArrayLike:
        """Estimates the shift that has occured between planned and telemetry sinograms.

        Parameters
        ----------
        mode : {"brute-force", "first-open", "center-of-mass"}, optional
            The mode parameter determines how the leaf shift in a given projection is
            estimated. Default is "first-open".
            "brute-force"
                Leaf shift is estimated based on least deviation between telemetry and
                all possible integer leaf adaptations of planned sinogram projections.
            "first-open"
                Leaf shift is estimated based on position of first open leaf.
            "last-open"
                Leaf shift is estimated based on positition of last open leaf.
            "center-of-mass"
                Leaf shift is estimated by center of mass comparison.
        flot_open_threshold : float, optional
            The threshold for fractional leaf open time to consider a leaf as being open.
            Default value is 0.2.

        Returns
        -------
        npt.ArrayLike
            Array containing estimated leaf adaptations for each projection of a sinogram.

        Notes
        -----
        Calculating the leaf shift is straightforward where the series of leaf open
        times are in agreement, with an index shift. However, for cases where there are
        variations in fluence (e.g., adaptation for a fraction of a projection, or
        erroneous leaf behaviour), more advanced registration approaches may be necessary.
        Where no leaves are open, it is not possible to tell if adaptation would have been
        correctly requested by the system.
        """
        if mode not in ["brute-force", "first-open", "last-open", "center-of-mass"]:
            raise ValueError("Invalid mode specified for leaf shift estimation.")
        self_data = self._df["channel_values"].to_numpy()
        other_data = other._df["channel_values"].to_numpy()
        if not self_data.shape == other_data.shape:
            raise ValueError(
                "Sinograms could not be compared with shapes",
                self_data.shape,
                other_data.shape,
            )
        estimated_shift = []
        for self_row, other_row in zip(self_data, other_data):
            if np.sum(self_row) < flot_open_threshold:
                estimated_shift.append(0)
            elif mode == "brute-force":
                pad_width = len(other_row)
                padded_target = np.pad(self_row, pad_width, mode="constant")
                best_shift = 0
                min_diff = float("inf")
                for shift in range(-pad_width, pad_width):
                    shifted_target = np.roll(padded_target, shift)
                    diff = np.sum(
                        (other_row - shifted_target[pad_width : 2 * pad_width]) ** 2
                    )
                    if diff < min_diff:
                        min_diff = diff
                        best_shift = shift
                estimated_shift.append(best_shift)
            elif mode == "first-open":
                estimated_shift.append(
                    np.argmax(other_row > flot_open_threshold)
                    - np.argmax(self_row > flot_open_threshold)
                )
            elif mode == "last-open":
                estimated_shift.append(
                    -(
                        (
                            len(self_row)
                            - 1
                            - np.argmax(self_row[::-1] > flot_open_threshold)
                        )
                        - (
                            len(other_row)
                            - 1
                            - np.argmax(other_row[::-1] > flot_open_threshold)
                        )
                    )
                )
            elif mode == "center-of-mass":
                estimated_shift.append(
                    round(center_of_mass(other_row)[0] - center_of_mass(self_row)[0])
                )
        return estimated_shift

    def detect_fluence_variations(
        self,
        other: RadixactSinogram,
        flot_threshold=0.1,
    ) -> npt.ArrayLike:
        """Detects projections containing adaptation (or errors).

        Parameters
        ----------
        other : RadixactSinogram
            The (telemetry or expected adaptation) sinogram against which to compare
            fractional leaf open times.
        flot_threshold : float, optional
            Threshold for fractional leaf open time variation. Default is 0.1.

        Returns
        -------
        npt.ArrayLike
            Sinogram row mask containg True and False values, where a fractional
            leaf open time exceeding the threshold occurs within a projection.
        """
        detected_adaptations = []
        for self_row, other_row in zip(
            self.fractional_leaf_open_times(),
            other.fractional_leaf_open_times(),
        ):
            if np.any(np.abs(self_row - other_row) > flot_threshold):
                detected_adaptations.append(True)
            else:
                detected_adaptations.append(False)
        return detected_adaptations

    def detect_fluence_errors(
        self,
        other: RadixactSinogram,
        error_detection_mode="integrated-leaf-open-time",
        total_flot_error_threshold=0.5,
    ) -> npt.ArrayLike:
        """Detects projections containing where fluence errors have occured.

        Parameters
        ----------
        other : RadixactSinogram
            The (telemetry) sinogram against which to compare leaf adaptations.
        error_detection_mode : {"integrated-leaf-open-time"}, optional
            Method to use to detect leaf adaptation errors. Default is
            "integrated-leaf-open-time".
        total_flot_error_threshold : float, optional
            Threshold for projection total fractional leaf open time. Default is 0.5.

        Returns
        -------
        npt.ArrayLike
            Sinogram row mask containing True and False values, for differences in
            total leaf open within a projection exceeding threshold.
        """
        detected_fluence_errors = []
        # TODO implement other approaches to this
        if error_detection_mode == "integrated-leaf-open-time":
            for self_row, other_row in zip(
                self._df["channel_values"].to_numpy(),
                other._df["channel_values"].to_numpy(),
            ):
                if (
                    abs(np.sum(self_row) - np.sum(other_row))
                    > total_flot_error_threshold
                ):
                    detected_fluence_errors.append(True)
                else:
                    detected_fluence_errors.append(False)
        return detected_fluence_errors

    def fractional_leaf_open_times(self) -> npt.ArrayLike:
        """Returns the sinogram fractional leaf open times.

        Returns
        -------
        ArrayLike
            Fractional leaf open times.
        """
        return self._df.get_column("channel_values").to_numpy()

    def leaf_open_times(self) -> npt.ArrayLike:
        """Returns the sinogram leaf open times.

        Returns
        -------
        ArrayLike
            Leaf open times in milliseconds.
        """
        projection_time = self.projection_duration()
        return self._df.get_column("channel_values").to_numpy() * projection_time

    def projection_duration(self) -> float:
        """Calculates the projection duration in milliseconds.

        Returns
        -------
        float
            Projection duration, in milliseconds.
        """
        return (
            self._df.get_column("time").to_numpy()[1]
            - self._df.get_column("time").to_numpy()[0]
        ) * 1000

    def to_csv(self, path: str | os.PathLike):
        """Writes sinogram leaf open times or channel signal data to CSV file.

        Parameters
        ----------
        path : str | os.PathLike
            Path of the CSV file to be written.
        """
        if self._df.shape[1] <= 64:
            header_row = ["Projection"] + [
                "Leaf %d LOT" % (i + 1) for i in range(self._df.shape[1])
            ]
        else:
            header_row = ["Projection"] + [
                "Channel %d count" % (i + 1) for i in range(self._df.shape[1])
            ]
        projection_column = np.arange(1, self._df.shape[0] + 1).reshape((-1, 1))
        np.savetxt(
            path,
            np.concatenate([projection_column, self._df], axis=1),
            fmt="%i" + (",%8f" * self._df.shape[1]),
            delimiter=",",
            header=",".join(header_row),
        )

    def to_dcm(self, path: str | os.PathLike, template_path: str | os.PathLike) -> None:
        """Writes sinogram leaf open times to DICOM RTPLAN file, using an existing
        DICOM RTPLAN file as a template.

        Parameters
        ----------
        path : str | os.PathLike
            Path of the DICOM RTPLAN file to be written.
        template_path : str | os.PathLike
            Path of the DICOM RTPLAN file to be used as a template.

        Notes
        -----
        Depending on how the DICOM RTPLAN has been collected, e.g., either exported
        from the TPS or collected from a PDE or delivery analysis cache, the first
        control point may or may not have associated leaf open times.
        """
        ds = pydicom.dcmread(template_path)
        if len(self._df) == len(ds.BeamSequence[0].ControlPointSequence) - 1:
            for cp, lot in zip(ds.BeamSequence[0].ControlPointSequence[0:-1], self._df):
                cp[0x300D, 0x10A7].value = "//".join([str(x) for x in lot]).encode(
                    "utf-8"
                )
        elif len(self._df) == len(ds.BeamSequence[0].ControlPointSequence) - 2:
            for cp, lot in zip(ds.BeamSequence[0].ControlPointSequence[1:-1], self._df):
                cp[0x300D, 0x10A7].value = "//".join([str(x) for x in lot]).encode(
                    "utf-8"
                )
        else:
            raise RuntimeError(
                "Template DICOM RTPLAN has incorrect number of control points"
            )
        ds.save_as(path)

    def to_parquet(self, path: str | os.PathLike) -> None:
        """Writes sinogram dataframe to parquet file.

        Parameters
        ----------
        path : str | os.PathLike
            Path of the parquet file to be written.
        """
        self._df.write_parquet(path)

    def upsample(self, frequency: int = 20) -> RadixactSinogram:
        """Upsamples sinogram to fractional projections.

        Parameters
        ----------
        frequency : int, optional
            Frequency of sampling, or number of discrete points per projection.
            Default is 20.

        Returns
        -------
        RadixactSinogram
            The sinogram upsampled to specified number of points per projection.

        Notes
        -----
        The function is intended to facilitate the prediction of leaf motion
        adaptations at a sub-projection interval. Upsampling of a previously
        upsampled sinogram is not currently validated. Upsampling of sinogram
        other than planned sinogram is not currently validated.
        """
        # Check if data is already downsampled, otherwise perform downsampling.
        curr_tau_interval = self._df["tau"][1] - self._df["tau"][0]
        if curr_tau_interval < 1:
            raise NotImplementedError(
                "Upsampling a previously upsampled sinogram is not validated"
            )
        # Calculate new projection, tau, time and gantry angles.
        curr_tau_start = self._df["tau"][0]
        curr_tau_end = self._df["tau"][-1] + curr_tau_interval
        curr_time_interval = self._df["time"][1] - self._df["time"][0]
        new_tau_interval = 1 / frequency
        new_tau = np.arange(curr_tau_start, curr_tau_end, new_tau_interval)
        new_projection = [int(x) for x in (new_tau - new_tau[0] + 1)]
        new_time = (new_tau - new_tau[0]) * curr_time_interval
        # new_time_interval = curr_time_interval / frequency
        # new_time = np.arange(curr_time_start, curr_time_end, new_time_interval)
        new_gantry_angles = ((new_tau % 51) * (360 / 51) - (180 / 51)) % 360
        # Create new data frame for updated leaf open times.
        data = pl.DataFrame(
            {
                "projection": new_projection,
                "tau": new_tau,
                "time": new_time,
                "gantry_angle": new_gantry_angles,
            }
        )
        # Calculate new leaf open times for resampled projections
        new_leaf_open_times = []
        for upsampled_row in data.iter_rows(named=True):
            upsampled_row_leaf_open_times = []
            downsampled_row = self._df.filter(
                pl.col("projection") == upsampled_row["projection"]
            )[0]
            upsampled_row_tau_start = upsampled_row["tau"] % 1
            upsampled_row_tau_end = (upsampled_row["tau"] + new_tau_interval) % 1
            upsampled_row_tau_interval = upsampled_row_tau_end - upsampled_row_tau_start
            # Calculate new leaf open time for leaf channel
            for downsampled_leaf_open_time in downsampled_row["channel_values"][0]:
                if downsampled_leaf_open_time == 0.0:
                    upsampled_row_leaf_open_times.append(0)
                else:
                    overlap_start = max(
                        0.5 - (downsampled_leaf_open_time / 2), upsampled_row_tau_start
                    )
                    overlap_end = min(
                        0.5 + (downsampled_leaf_open_time / 2), upsampled_row_tau_end
                    )
                    overlap = (
                        overlap_end - overlap_start
                        if overlap_start <= overlap_end
                        else 0
                    )
                    upsampled_row_leaf_open_times.append(
                        overlap / upsampled_row_tau_interval
                    )
            new_leaf_open_times.append(upsampled_row_leaf_open_times)
        data = data.with_columns(
            channel_values=pl.Series(np.array(new_leaf_open_times))
        )
        return type(self)(data)

    # endregion

    # region Private methods

    # endregion

    # region Static methods

    @staticmethod
    def plot_leaf_adaptation(
        plan_sinogram: RadixactSinogram,
        telemetry_sinogram: RadixactSinogram,
        telemetry_timing: RadixactTiming,
        motion: RadixactSynchronyMotion,
        projection_limit=None,
        leaf_limit=None,
        leaf_adaptation_threshold: float = 4.3 / 6.25,
        error_threshold: float = 0.5,
        registration_mode: str = "brute-force",
        planned_sinogram_colormap: str = "Greys",
        expected_sinogram_colormap: str = "Blues",
        telemetry_sinogram_colormap: str = "Greens",
        erroneous_sinogram_colormap: str = "Oranges",
        observed_adaptation_label: str = "Telemetry adaptation",
        closed_leaf_threshold: float = 0.05,
        closed_leaf_label: str = "Closed leaves",
        closed_leaf_color: str = "lightgrey",
        closed_leaf_alpha: float = 0.3,
        interrupt_label: str = "Beam interrupt",
        interrupt_color: str = "brown",
        interrupt_shown: bool = True,
        figsize=(12, 14),
    ):
        # TODO: assert agreement in # projections
        projections = len(plan_sinogram.fractional_leaf_open_times())
        if projection_limit is None:
            proj_limit_lower = 0
            proj_limit_upper = projections
        else:
            proj_limit_lower = projection_limit[0]
            proj_limit_upper = projection_limit[1]
        if leaf_limit is None:
            leaf_limit_lower = 1
            leaf_limit_upper = 64
        else:
            leaf_limit_lower = leaf_limit[0]
            leaf_limit_upper = leaf_limit[1]
        if (leaf_limit_upper - leaf_limit_lower) > 36:
            multiple_locator_divisor = 8
        else:
            multiple_locator_divisor = 4
        # Calculate expected sinogram
        optimal_sinogram = plan_sinogram.apply_motion_adaptation(
            telemetry_timing, motion, leaf_adaptation_mode="optimal"
        )
        upsampled_sinogram = plan_sinogram.upsample().apply_motion_adaptation(
            telemetry_timing,
            motion,
            leaf_adaptation_threshold=leaf_adaptation_threshold,
        )
        expected_sinogram = upsampled_sinogram.downsample()
        # Calculate leaf shift
        estimated_leaf_shifts = plan_sinogram.estimated_leaf_adaptation(
            telemetry_sinogram, mode=registration_mode
        )
        # Identify closed leaves
        closed_leaves = plan_sinogram.closed_leaf_projections(closed_leaf_threshold)
        # Detect fluence variation
        variation = plan_sinogram.detect_fluence_errors(
            telemetry_sinogram, total_flot_error_threshold=error_threshold
        )
        # Calculate pause time
        interrupts = np.unique(
            (
                np.array(telemetry_timing._df["tau"][1:])
                - telemetry_timing._df["tau"][0]
            )[:-1]
        )
        # Get image data
        planned_im_data = plan_sinogram.fractional_leaf_open_times()[
            proj_limit_lower:proj_limit_upper
        ]
        # detector_im_data = detector_sinogram.image_data()[proj_limit_lower:proj_limit_upper]
        # detector_im_data = detector_sinogram.image_data()
        telemetry_im_data = telemetry_sinogram.fractional_leaf_open_times()[
            proj_limit_lower:proj_limit_upper
        ]
        expected_im_data = expected_sinogram.fractional_leaf_open_times()[
            proj_limit_lower:proj_limit_upper
        ]
        # Define plot grid
        fig, axs = plt.subplots(nrows=6, ncols=1, sharex=True, figsize=figsize)

        # Define paired x axis functions
        def proj_to_rotation(x):
            return x / 51

        def rotation_to_proj(x):
            return x * 51

        # Plot planned sinogram
        planned_im = axs[0].imshow(
            np.rot90(planned_im_data)[leaf_limit_lower - 1 : leaf_limit_upper],
            extent=(
                proj_limit_lower,
                proj_limit_upper,
                leaf_limit_lower,
                leaf_limit_upper,
            ),
            vmin=0,
            vmax=1,
            aspect="auto",
            cmap=mpl.colormaps[planned_sinogram_colormap],
            interpolation="none",
        )
        axs[0].set_ylabel("Leaf")
        axs[0].yaxis.set_major_locator(plt.MultipleLocator(multiple_locator_divisor))
        secax0 = axs[0].secondary_xaxis(
            "top", functions=(proj_to_rotation, rotation_to_proj)
        )
        secax0.xaxis.set_major_locator(ticker.MultipleLocator(1))
        secax0.set_xlabel("Gantry rotations")
        axs[0].legend(
            loc="lower center",
            handles=[
                mpatches.Patch(
                    color=mpl.colormaps[planned_sinogram_colormap](0.8),
                    label="Planned sinogram",
                )
            ],
        )
        planned_cax = axs[0].inset_axes([1.01, 0, 0.01, 1])
        planned_cb = fig.colorbar(planned_im, cax=planned_cax, format="%.1f")
        planned_cb.set_label("Leaf open time")
        # Plot target offset
        axs[1].plot(
            optimal_sinogram._df["projection"],
            optimal_sinogram._df["target_offset_beams_eye_view"],
            color="C1",
            label="Beam's eye view",
            linewidth=1.5,
            zorder=4,
        )
        axs[1].plot(
            optimal_sinogram._df["projection"],
            optimal_sinogram._df["target_offset_x"],
            color="C7",
            label="IEC-X (left-right)",
            linewidth=1,
            zorder=3,
        )
        axs[1].plot(
            optimal_sinogram._df["projection"],
            optimal_sinogram._df["target_offset_y"],
            color="C9",
            label="IEC-Y (superior-inferior)",
            linewidth=1,
            zorder=2,
        )
        axs[1].plot(
            optimal_sinogram._df["projection"],
            optimal_sinogram._df["target_offset_z"],
            color="C8",
            label="IEC-Z (anterior-posterior)",
            linewidth=1,
            zorder=1,
        )
        axs[1].legend(loc="lower center", ncol=4)
        axs[1].sharex(axs[0])
        axs[1].set_ylabel("Target offset (mm)")
        secax1 = axs[1].secondary_xaxis(
            "top", functions=(proj_to_rotation, rotation_to_proj)
        )
        secax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
        secax1.xaxis.set_major_formatter(ticker.NullFormatter())
        offset_ylim = (
            round(
                max(
                    max(np.abs(optimal_sinogram._df["target_offset_x"])),
                    max(np.abs(optimal_sinogram._df["target_offset_y"])),
                    max(np.abs(optimal_sinogram._df["target_offset_z"])),
                )
                / 6.25
            )
            + 1
        )
        axs[1].set_ylim(-6.25 * offset_ylim, 6.25 * offset_ylim)
        # Plot adaptation
        axs[2].step(
            optimal_sinogram._df["projection"],
            optimal_sinogram._df["leaf_shift"],
            color="C1",
            label="Beam's eye view",
            linewidth=1.5,
            where="post",
            zorder=3,
        )
        axs[2].step(
            upsampled_sinogram._df["projection"],
            upsampled_sinogram._df["leaf_shift"],
            color=mpl.colormaps[expected_sinogram_colormap](0.8),
            label="Expected adaptation",
            linewidth=1.5,
            where="mid",
            zorder=2,
        )
        axs[2].fill_between(
            plan_sinogram._df["projection"],
            estimated_leaf_shifts,
            color=mpl.colormaps[telemetry_sinogram_colormap](0.3),
            step="post",
            label=observed_adaptation_label,
            zorder=1,
        )
        adaptation_ylim = (
            round(max(np.abs(optimal_sinogram._df["leaf_shift"].to_list()))) + 1
        )
        if True in closed_leaves:
            axs[2].fill_between(
                plan_sinogram._df["projection"],
                [-adaptation_ylim] * projections,
                [adaptation_ylim] * projections,
                where=closed_leaves,
                step="post",
                alpha=closed_leaf_alpha,
                color=closed_leaf_color,
                label=closed_leaf_label,
            )
        adaptation_ticks = range(-adaptation_ylim, adaptation_ylim + 1)
        if len(interrupts) > 0:
            adaptation_ncol = 4
        else:
            adaptation_ncol = 3
        axs[2].legend(loc="lower center", ncol=adaptation_ncol)
        axs[2].set_xlim(proj_limit_lower, proj_limit_upper)
        axs[2].set_ylabel("Leaf adaptation")
        axs[2].set_ylim(-adaptation_ylim, adaptation_ylim)
        axs[2].set_yticks(
            adaptation_ticks, labels=[str(tick) for tick in adaptation_ticks]
        )
        secax2 = axs[2].secondary_xaxis(
            "top", functions=(proj_to_rotation, rotation_to_proj)
        )
        secax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
        secax2.xaxis.set_major_formatter(ticker.NullFormatter())
        # Plot expected sinogram
        axs[3].imshow(
            np.rot90(planned_im_data)[leaf_limit_lower - 1 : leaf_limit_upper],
            extent=(
                proj_limit_lower,
                proj_limit_upper,
                leaf_limit_lower,
                leaf_limit_upper,
            ),
            vmin=0,
            vmax=1,
            aspect="auto",
            cmap=mpl.colormaps[planned_sinogram_colormap],
            interpolation="none",
        )
        expected_row_mask = (
            np.abs(planned_im_data - expected_im_data) < error_threshold
        ).all(axis=1)
        expected_full_mask = np.tile(
            expected_row_mask[:, np.newaxis], (1, expected_im_data.shape[1])
        )
        expected_masked_sinogram = np.ma.masked_array(
            expected_im_data, mask=expected_full_mask
        )
        expected_overlay_im = axs[3].imshow(
            np.rot90(expected_masked_sinogram)[leaf_limit_lower - 1 : leaf_limit_upper],
            extent=(
                proj_limit_lower,
                proj_limit_upper,
                leaf_limit_lower,
                leaf_limit_upper,
            ),
            vmin=0,
            vmax=1,
            aspect="auto",
            cmap=mpl.colormaps[expected_sinogram_colormap],
            interpolation="none",
        )
        expected_handles = [
            mpatches.Patch(
                color=mpl.colormaps[planned_sinogram_colormap](0.8),
                label="Planned sinogram",
            ),
            mpatches.Patch(
                color=mpl.colormaps[expected_sinogram_colormap](0.8),
                label="Expected sinogram variation",
            ),
        ]
        expected_labels = ["Planned sinogram", "Expected sinogram variation"]
        axs[3].set_ylabel("Leaf")
        axs[3].yaxis.set_major_locator(plt.MultipleLocator(multiple_locator_divisor))
        axs[3].legend(
            loc="lower center", handles=expected_handles, labels=expected_labels, ncol=4
        )
        expected_variations = 1 - np.sum(expected_row_mask) / len(expected_row_mask)
        axs[3].text(
            1,
            0.9,
            f"Expected projections with variations = {expected_variations * 100:.0f}% ",
            horizontalalignment="right",
            verticalalignment="center",
            transform=axs[3].transAxes,
        )
        secax3 = axs[3].secondary_xaxis(
            "top", functions=(proj_to_rotation, rotation_to_proj)
        )
        secax3.xaxis.set_major_locator(ticker.MultipleLocator(1))
        secax3.xaxis.set_major_formatter(ticker.NullFormatter())
        expected_cax = axs[3].inset_axes([1.01, 0, 0.01, 1])
        expected_cb = fig.colorbar(expected_overlay_im, cax=expected_cax, format="%.1f")
        expected_cb.set_label("Leaf open time")
        # Plot telemetry sinogram.
        axs[4].imshow(
            np.rot90(planned_im_data)[leaf_limit_lower - 1 : leaf_limit_upper],
            extent=(
                proj_limit_lower,
                proj_limit_upper,
                leaf_limit_lower,
                leaf_limit_upper,
            ),
            vmin=0,
            vmax=1,
            aspect="auto",
            cmap=mpl.colormaps[planned_sinogram_colormap],
            interpolation="none",
        )
        telemetry_row_mask = (
            np.abs(planned_im_data - telemetry_im_data) < error_threshold
        ).all(axis=1)
        telemetry_full_mask = np.tile(
            telemetry_row_mask[:, np.newaxis], (1, telemetry_im_data.shape[1])
        )
        telemetry_masked_sinogram = np.ma.masked_array(
            telemetry_im_data, mask=telemetry_full_mask
        )
        telemetry_overlay_im = axs[4].imshow(
            np.rot90(telemetry_masked_sinogram)[
                leaf_limit_lower - 1 : leaf_limit_upper
            ],
            extent=(
                proj_limit_lower,
                proj_limit_upper,
                leaf_limit_lower,
                leaf_limit_upper,
            ),
            vmin=0,
            vmax=1,
            aspect="auto",
            cmap=mpl.colormaps[telemetry_sinogram_colormap],
            interpolation="none",
        )
        handles = [
            mpatches.Patch(
                color=mpl.colormaps[planned_sinogram_colormap](0.8),
                label="Planned sinogram",
            ),
            mpatches.Patch(
                color=mpl.colormaps[telemetry_sinogram_colormap](0.8),
                label="Telemetry sinogram variation",
            ),
        ]
        labels = ["Planned sinogram", "Telemetry sinogram variation"]
        erroneous_row_mask = (
            np.abs(np.sum(planned_im_data, axis=1) - np.sum(telemetry_im_data, axis=1))
            < error_threshold
        )
        erroneous_full_mask = np.tile(
            erroneous_row_mask[:, np.newaxis], (1, telemetry_im_data.shape[1])
        )
        erroneous_masked_sinogram = np.ma.masked_array(
            telemetry_im_data, mask=erroneous_full_mask
        )
        axs[4].imshow(
            np.rot90(erroneous_masked_sinogram)[
                leaf_limit_lower - 1 : leaf_limit_upper
            ],
            extent=(
                proj_limit_lower,
                proj_limit_upper,
                leaf_limit_lower,
                leaf_limit_upper,
            ),
            vmin=0,
            vmax=1,
            aspect="auto",
            cmap=mpl.colormaps[erroneous_sinogram_colormap],
            interpolation="none",
        )
        if not erroneous_row_mask.all():
            handles.append(
                mpatches.Patch(
                    color=mpl.colormaps[erroneous_sinogram_colormap](0.8),
                    label="Erroneous sinogram variation",
                )
            )
            labels.append("Erroneous sinogram variation")
        if interrupt_shown:
            for interrupt in interrupts:
                axs[4].vlines(
                    x=interrupt,
                    ymin=leaf_limit_lower,
                    ymax=leaf_limit_upper,
                    color=interrupt_color,
                    label=interrupt_label,
                )
                handle, label = axs[4].get_legend_handles_labels()
                handles.append(handle[0])
                labels.append(label[0])
        axs[4].set_ylabel("Leaf")
        axs[4].yaxis.set_major_locator(plt.MultipleLocator(multiple_locator_divisor))
        axs[4].legend(loc="lower center", handles=handles, labels=labels, ncol=3)
        telemetry_variations = 1 - np.sum(telemetry_row_mask) / len(telemetry_row_mask)
        axs[4].text(
            1,
            0.9,
            f"Projections with variations = {telemetry_variations * 100:.0f}% ",
            horizontalalignment="right",
            verticalalignment="center",
            transform=axs[4].transAxes,
        )
        secax4 = axs[4].secondary_xaxis(
            "top", functions=(proj_to_rotation, rotation_to_proj)
        )
        secax4.xaxis.set_major_locator(ticker.MultipleLocator(1))
        secax4.xaxis.set_major_formatter(ticker.NullFormatter())
        telemetry_cax = axs[4].inset_axes([1.01, 0, 0.01, 1])
        telemetry_cb = fig.colorbar(
            telemetry_overlay_im, cax=telemetry_cax, format="%.1f"
        )
        telemetry_cb.set_label("Leaf open time")
        # Plot total leaf open time
        axs[5].step(
            plan_sinogram._df["projection"],
            np.sum(plan_sinogram.fractional_leaf_open_times(), axis=1),
            color=mpl.colormaps[planned_sinogram_colormap](1.0),
            label="Planned sinogram",
            linewidth=1.5,
            where="mid",
        )
        axs[5].step(
            plan_sinogram._df["projection"],
            np.sum(telemetry_sinogram.fractional_leaf_open_times(), axis=1),
            color=mpl.colormaps[telemetry_sinogram_colormap](0.6),
            label="Telemetry sinogram",
            linewidth=1.5,
            where="mid",
        )
        total_lot_ylim = (
            int(
                max(
                    np.max(np.sum(plan_sinogram.fractional_leaf_open_times(), axis=1)),
                    np.max(
                        np.sum(telemetry_sinogram.fractional_leaf_open_times(), axis=1)
                    ),
                )
            )
            + 1
        )
        if True in closed_leaves:
            axs[5].fill_between(
                plan_sinogram._df["projection"],
                0,
                total_lot_ylim,
                where=closed_leaves,
                step="post",
                alpha=closed_leaf_alpha,
                color=closed_leaf_color,
                label=closed_leaf_label,
            )
        if True in variation:
            axs[5].fill_between(
                plan_sinogram._df["projection"],
                0,
                total_lot_ylim,
                where=variation,
                step="post",
                color=mpl.colormaps[erroneous_sinogram_colormap](0.1),
                label="Erroneous variation",
            )
        axs[5].set_ylim(0, total_lot_ylim)
        axs[5].set_ylabel("Total leaf open time")
        axs[5].set_xlabel("Projection")
        secax5 = axs[5].secondary_xaxis(
            "top", functions=(proj_to_rotation, rotation_to_proj)
        )
        secax5.xaxis.set_major_locator(ticker.MultipleLocator(1))
        secax5.xaxis.set_major_formatter(ticker.NullFormatter())
        if interrupt_shown:
            for interrupt in interrupts:
                axs[5].vlines(
                    x=interrupt,
                    ymin=0,
                    ymax=total_lot_ylim,
                    color=interrupt_color,
                    label=interrupt_label,
                    linewidth=1.5,
                )
        axs[5].legend(loc="lower center", ncol=5)
        # Return figure
        fig.tight_layout()
        return fig

    # end region
