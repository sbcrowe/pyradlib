# -*- coding: utf-8 -*-
"""Motion analysis module.

This module provides functionality for processing of motion data from Synchrony
treatments.
"""

from __future__ import annotations

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import logging
import os
import xml.etree.ElementTree as et
from functools import cached_property

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import numpy.typing as npt
import polars as pl

logger = logging.getLogger(__name__)


class RadixactSynchronyMotion:
    # region Constructors

    def __init__(self, df: pl.DataFrame, uid: str = None) -> RadixactSynchronyMotion:
        """
        Parameters
        ----------
        df : pl.DataFrame
            Dataframe containing relative timestamps of motion data (in sec),
            IEC-X, IEC-Y, and IEC-Z displacements (in mm), and optionally,
            rigid body values, potential diffs and measured deltas (in mm).
        uid : str, optional
            Unique identifier string, to allow association with other data. Default is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
        self._df = df
        self._uid = uid

    @classmethod
    def from_array_likes(
        cls,
        times: npt.ArrayLike,
        target_offset_x: npt.ArrayLike,
        target_offset_y: npt.ArrayLike,
        target_offset_z: npt.ArrayLike,
        rigid_bodys: npt.ArrayLike = None,
        potential_diffs: npt.ArrayLike = None,
        measured_diffs: npt.ArrayLike = None,
        uid: str = None,
    ) -> RadixactSynchronyMotion:
        """_summary_

        Parameters
        ----------
        times : npt.ArrayLike
            Relative timestamps of motion data points, in sec.
        target_offset_x : npt.ArrayLike
            IEC-X or left-right displacements, in mm.
        target_offset_y : npt.ArrayLike
            IEC-Y or superior-inferior displacements, in mm.
        target_offset_z : npt.ArrayLike
            IEC-Z or anterior-posterior displacements, in mm.
        rigid_bodys : npt.ArrayLike, optional
            Rigid body fiducial deformation, in mm.
        potential_diffs : npt.ArrayLike, optional
            Potential difference of uncertainty in model, in mm.
        measured_diffs : npt.ArrayLike, optional
            Measured difference between model and images, in mm.
        uid : str, optional
            Unique identifier string, to allow association with other data. Default is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
        # TODO add timestamp reading
        df_dict = {}
        nan_indices = np.where(np.isnan(times))[0]
        # TODO Verify segment_index and not result_index is appropriate below.
        if len(nan_indices) == 0:
            segment = [1] * len(times)
        else:
            segment = []
            curr_fragment = 1
            curr_index = 0
            for nan_index in nan_indices:
                num_fragment_data_points = nan_index - curr_index
                segment += [curr_fragment] * num_fragment_data_points
                segment += [np.nan]
                curr_index = nan_index + 1
                curr_fragment += 1
            segment.pop()
        df_dict["segment_index"] = segment
        df_dict["delta_time"] = times
        df_dict["target_offset_x"] = target_offset_x
        df_dict["target_offset_y"] = target_offset_y
        df_dict["target_offset_z"] = target_offset_z
        if rigid_bodys is not None:
            df_dict["rigid_body"] = rigid_bodys
        if potential_diffs is not None:
            df_dict["potential_diff"] = potential_diffs
        if measured_diffs is not None:
            df_dict["measured_diff"] = measured_diffs
        df = pl.DataFrame(df_dict)
        df = cls._calculate_additional_offsets(df)
        return cls(df, uid)

    @classmethod
    def from_motion_extractor_csv(
        cls,
        path_x: str | os.PathLike,
        path_y: str | os.PathLike,
        path_z: str | os.PathLike,
        uid: str = None,
    ) -> RadixactSynchronyMotion:
        """Extracts motion data from csv file representing a patient
        treatment fraction fragment, taken from Motion Data Extractor tool.

        Parameters
        ----------
        path_x : str | os.PathLike
            Path to csv file containing IEC-X or left-right motion data.
        path_y : str | os.PathLike
            Path to csv file containing IEC-Y or superior-inferior motion data.
        path_z : str | os.PathLike
            Path to csv file containing IEC-Z or anterior-posterior motion data.
        uid : str, optional
            Unique identifier string, to allow association with other data. Defualt is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony Motion data encapsulated with helper functions.
        """
        times = np.loadtxt(path_x, delimiter=",", skiprows=1)[:, 0]
        x_displacements = np.loadtxt(path_x, delimiter=",", skiprows=1)[:, 1]
        y_displacements = np.loadtxt(path_y, delimiter=",", skiprows=1)[:, 1]
        z_displacements = np.loadtxt(path_z, delimiter=",", skiprows=1)[:, 1]
        return cls.from_array_likes(
            times, x_displacements, y_displacements, z_displacements, uid=uid
        )

    @classmethod
    def from_npz(
        cls, path: str | os.PathLike, uid: str = None
    ) -> RadixactSynchronyMotion:
        """Reads motion data stored in numpy npz file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the npz file containing motion data.
        uid : str, optional
            Unique identifier string, to allow association with other data. Default is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
        with np.load(path) as npz_data:
            df = pl.DataFrame(dict(npz_data))
            df = cls._calculate_additional_offsets(df)
            df = cls._calculate_additional_times(df)
            return cls(df, uid)

    @classmethod
    def from_parquet(
        cls, path: str | os.PathLike, uid: str = None
    ) -> RadixactSynchronyMotion:
        """Reads motion data stored in parquet file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to parquet file containing motion data.
        uid : str, optional
            Unique identifier string, to allow association with other data. Default is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
        df = pl.read_parquet(path)
        return cls(df, uid)

    @classmethod
    def from_patient_motions(
        cls, motions: list[RadixactSynchronyMotion]
    ) -> RadixactSynchronyMotion:
        """Combines motion data from multiple patients, typically all motions for
        a patient dataset cohort.

        Parameters
        ----------
        motions : list[RadixactSynchronyMotion]
            The patient motions to be combined.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
        result = pl.concat(
            [
                motion._df.with_columns(patient_index=pl.lit(key))
                for key, motion in enumerate(motions)
            ]
        ).select([pl.col("patient_index"), pl.all().exclude("patient_index")])
        logger.info(f"Combined {len(motions)} motion data files")
        return RadixactSynchronyMotion(result)

    @classmethod
    def from_session_motions(
        cls, motions: list[RadixactSynchronyMotion]
    ) -> RadixactSynchronyMotion:
        """Combines motion data from multiple sessions, typically all motions for
        a patient dataset.

        Parameters
        ----------
        motions : list[RadixactSynchronyMotion]
            The session motions to be combined.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
        result = pl.concat(
            [
                motion._df.with_columns(session_index=pl.lit(key))
                for key, motion in enumerate(motions)
            ]
        ).select([pl.col("session_index"), pl.all().exclude("session_index")])
        logger.info(f"Combined {len(motions)} motion data files")
        return RadixactSynchronyMotion(result)

    @classmethod
    def from_xml(
        cls, path: str | os.PathLike, uid: str = None
    ) -> RadixactSynchronyMotion:
        """Extracts motion data from *motionData.xml file representing a patient
        treatment fraction fragment, taken from a Delivery Analysis cache or Patient
        Data Extractor archive.

        Parameters
        ----------
        path : str | os.PathLike
            Path to XML file containing motion data.
        uid : str, optional
            Unique identifier string, to allow association with other data. Default is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.

        Notes
        -----
        Motion data are included in Patient Data Extractor archives, though are renamed
        using tokens described within the PatientExportDataBO.xml file. These files are
        also cached as C:/tomo/da/pts/URnumber/*motionData.xml when patient data is
        loaded within the Delivery Analysis tool. There may be multiple files per
        fraction.

        This function has only been validated for version 3 MotionData.xml files.

        Pauses in fractions (e.g. due to user, Synchrony, or system errors) will be
        indicated by numpy.nan. This allows more accurate plotting of the data using
        matplotlib.
        """
        tree = et.parse(path)
        root = tree.getroot()
        dfs = []
        for delivery_segment in root.iter("DeliverySegment"):
            for radiation_results in delivery_segment.iter("RadiationResults"):
                for model_data in radiation_results.iter("ModelData"):
                    data = []
                    for datapoint in model_data.iter("DataPoint"):
                        data.append(
                            {
                                "segment_index": int(delivery_segment[0].text),
                                "result_index": int(radiation_results[0].text),
                                "potential_diff_tolerance": float(model_data[0].text),
                                "measured_diff_tolerance": float(model_data[1].text),
                                "timestamp": int(datapoint[0].text),
                                "potential_diff": float(datapoint[1].text),
                                "rigid_body": float(datapoint[2].text),
                                "target_offset_x": float(datapoint[3][0].text),
                                "target_offset_y": float(datapoint[3][1].text),
                                "target_offset_z": float(datapoint[3][2].text),
                                "measured_diff": float(datapoint[4].text),
                            }
                        )
                    if len(data) > 0:
                        df = pl.DataFrame(data)
                        df = cls._calculate_additional_offsets(df)
                        dfs.append(df)
        if len(dfs) == 0:
            return None
        # Create dataframe
        df = pl.concat(dfs)
        df = cls._calculate_additional_times(df)
        # Remove any non-existent data
        df = df.fill_nan(None)
        df_cleaned = df.select(
            [pl.col(col) for col in df.columns if not df[col].is_null().all()]
        )
        return cls(df_cleaned, uid)

    # endregion

    # region Magic methods

    def __add__(self, other: RadixactSynchronyMotion) -> RadixactSynchronyMotion:
        """Concatenates two motion datasets together.

        Parameters
        ----------
        other : RadixactSynchronyMotion
            The Radixact Synchrony Motion data to be concated to this motion data.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony Motion data encapsulated with helper functions.

        Raises
        ------
        TypeError
            If either self or other is not Radixact Synchrony motion data.
        """
        if not isinstance(other, type(self)):
            raise TypeError("Unsupported operand type for +")
        return type(self)(pl.concat([self._df, other._df]))

    def __len__(self) -> int:
        """Calculates the number of points in the motion data as length.

        Returns
        -------
        int
            Length of motion data, or number of discrete points.
        """
        return len(self._df)

    def __radd__(self, other: RadixactSynchronyMotion) -> RadixactSynchronyMotion:
        """Concatenates two motion datasets together, for sum() functions.

        Parameters
        ----------
        other : RadixactSynchronyMotion
            The Radixact Synchrony Motion data to be concatenated to this motion data.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony Motion data encapsulated with helper functions.
        """
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactSynchronyMotion(uid={self._uid})"

    # endregion

    # region Properties

    @cached_property
    def metrics(self) -> pl.DataFrame:
        """Calculates motion and adaptation metrics.

        Returns
        -------
        pl.DataFrame
            The calculated motion and adaptation metrics.
        """
        if "patient_index" in self._df:
            result_df = []
            for patient_index, patient_df in self._df.group_by("patient_index"):
                for session_index, session_df in patient_df.group_by("session_index"):
                    result_df.append(
                        self._calculate_metrics(session_df).with_columns(
                            pl.lit(patient_index[0]).alias("patient_index"),
                            pl.lit(session_index[0]).alias("session_index"),
                        )
                    )
                result_df.append(
                    self._calculate_metrics(session_df).with_columns(
                        pl.lit(patient_index[0]).alias("patient_index"),
                        pl.lit(None).alias("session_index"),
                    )
                )
            result_df.append(
                self._calculate_metrics(session_df).with_columns(
                    pl.lit(None).alias("patient_index"),
                    pl.lit(None).alias("session_index"),
                )
            )
            return pl.concat(result_df).select(
                [
                    pl.col("patient_index"),
                    pl.col("session_index"),
                    pl.all().exclude(["patient_index", "session_index"]),
                ]
            )
        elif "session_index" in self._df:
            result_df = []
            for session_index, session_df in self._df.group_by("session_index"):
                result_df.append(
                    self._calculate_metrics(session_df).with_columns(
                        pl.lit(session_index[0]).alias("session_index")
                    )
                )
            result_df.append(
                self._calculate_metrics(session_df).with_columns(
                    pl.lit(None).alias("session_index"),
                )
            )
            return pl.concat(result_df).select(
                [pl.col("session_index"), pl.all().exclude("session_index")]
            )
        else:
            return self._calculate_metrics(self._df)

    # endregion

    # region Public methods

    def patient_fraction_less_than_threshold(
        self, offset_type: str = "target_offset_vector", threshold_step: float = 1
    ) -> pl.DataFrame:
        """Calculates the fraction of patient motion statistics less than thresholds,
        intended to provide insight on margin selection.

        Parameters
        ----------
        offset_type : str, optional
            The motion metric percentiles to evaluate against thresholds. Default is
            "target_offset_vector".
        threshold_step : float, optional
            Step width to use to define thresholds against which to evaluate offset
            percentiles, in mm. Default is 1.

        Returns
        -------
        pl.DataFrame
            DataFrame containing defining fraction of patient offset percentiles less
            than calculated thresholds.

        Notes
        -----
        This calculation is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        metrics = self.metrics.filter(
            pl.col("session_index").is_null(), pl.col("patient_index").is_not_null()
        )
        return self._fraction_less_than_threshold(metrics, offset_type, threshold_step)

    def plot_motion_histogram(
        self,
        mode: str = "absolute",
        fig_size: tuple[float, float] = (12, 4),
        offset_lim: tuple[float, float] = None,
        offset_bin: float = 0.25,
        vector_lim: tuple[float, float] = (0, 10),
        vector_bin: float = 0.25,
        title: str = None,
    ) -> mpl.Figure:
        """Generate a histogram plot of distribution of offset values in each dimension.

        Parameters
        ----------
        mode : str, optional
            Indicates whether to use absolute or relative offset. Default is "absolute".
        fig_size : tuple[float, float], optional
            Size of the figure. Default is (12, 4) inches.
        offset_lim : tuple[float, float], optional
            Minimum and maximum bin boundaries for offset. Default is None, in which
            case the limits will be calculated with floor(min) and ceiling(max),
            respectively.
        offset_bin : float, optional
            Bin width. Default is 0.25 mm.
        vector_lim : tuple[float, float], optional
            Minimum and maximum bin boundaries for vector. Default is None, in which
            case the limits will be 0 and ceiling(max).
        vector_bin : float, optional
            Bin width, Default is 0.25 mm.
        title : str, optional
            Title for figure. Default is None.

        Returns
        -------
        mpl.Figure
            Histogram of target offset values in each dimension.
        """
        if mode == "absolute":
            offsets = self._df.select(
                "target_offset_x",
                "target_offset_y",
                "target_offset_z",
                "target_offset_vector",
            ).to_numpy()
        if mode == "delta" or mode == "relative":
            offsets = self._df.select(
                "delta_target_offset_x",
                "delta_target_offset_y",
                "delta_target_offset_z",
                "delta_target_offset_vector",
            ).to_numpy()
        fig, axs = plt.subplots(
            ncols=4,
            nrows=1,
            gridspec_kw={"hspace": 0, "wspace": 0},
            constrained_layout=True,
            figsize=fig_size,
        )
        if offset_lim is None:
            min_offset = np.floor(np.min(offsets[:, 0:3]))
            max_offset = np.ceil(np.max(offsets[:, 0:3]))
        else:
            min_offset = offset_lim[0]
            max_offset = offset_lim[1]
        for index, axis in enumerate(["x", "y", "z"]):
            # TODO Don't assume that max_offset - min_offset is divisible by offset_bin
            y, x, _ = axs[index].hist(
                np.clip(offsets[:, index], min_offset, max_offset),
                list(
                    np.arange(
                        min_offset,
                        max_offset + offset_bin,
                        offset_bin,
                    )
                ),
                density=False,
            )
            axs[index].axvline(
                np.nanmean(offsets[:, index]),
                color="r",
                linestyle="-",
                label=r"$\bar{"
                + axis
                + r"}$"
                + f"={np.nanmean(offsets[:, index]):0.1f}\n"
                + r"$2\sigma_{"
                + axis
                + r"}$"
                + f"={2 * np.nanstd(offsets[:, index]):0.1f}",
            )
            axs[index].axvline(
                np.nanmean(offsets[:, index]) - 2 * np.nanstd(offsets[:, index]),
                color="r",
                linestyle="--",
                label="_",
            )
            axs[index].axvline(
                np.nanmean(offsets[:, index]) + 2 * np.nanstd(offsets[:, index]),
                color="r",
                linestyle="--",
                label="_",
            )
            axs[index].legend(
                loc="upper right", handlelength=0, handletextpad=0, labelcolor="r"
            )
            axs[index].set_xlim((min_offset, max_offset))
            axs[index].ticklabel_format(scilimits=[-3, 3])
        if mode == "absolute":
            axs[0].set_title("IEC-X (left-right)")
            axs[1].set_title("IEC-Y (sup-inf)")
            axs[2].set_title("IEC-Z (ant-post)")
        elif mode == "delta" or mode == "relative":
            axs[0].set_title("ΔIEC-X (left-right)")
            axs[1].set_title("ΔIEC-Y (sup-inf)")
            axs[2].set_title("ΔIEC-Z (ant-post)")
        if vector_lim is None:
            min_vector = 0
            max_vector = np.ceil(np.max(offsets[:, 3]))
        else:
            min_vector = vector_lim[0]
            max_vector = vector_lim[1]
        axs[3].hist(
            np.clip(offsets[:, 3], min_vector, max_vector),
            list(
                np.arange(
                    min_vector,
                    max_vector + vector_bin,
                    vector_bin,
                )
            ),
            density=False,
            cumulative=-1,
        )
        axs[3].axvline(
            np.nanpercentile(offsets[:, 3], 95),
            color="g",
            linestyle="-",
            label="$r_{95}$" + f"={np.percentile(offsets[:, 3], 95):0.1f}",
        )
        axs[3].legend(
            loc="upper right", handlelength=0, handletextpad=0, labelcolor="g"
        )
        axs[3].set_xlim((min_vector, max_vector))
        axs[3].ticklabel_format(scilimits=[-3, 3])
        if mode == "absolute":
            axs[3].set_title("r (vector)")
        elif mode == "delta" or mode == "relative":
            axs[3].set_title("Δr (vector)")
        fig.supylabel("Number of data points")
        if mode == "absolute":
            fig.supxlabel("Target offset relative to position at registration (mm)")
        elif mode == "delta" or mode == "relative":
            fig.supxlabel(
                "Target offset relative to position at start of tracking (mm)"
            )
        if title is not None:
            fig.suptitle(title)
        return fig

    def plot_patient_fraction_less_than_threshold(
        self,
        offset_type: str = "target_offset_vector",
        threshold_step: float = 1,
        figsize=(12, 4),
    ) -> mpl.Figure:
        """Plots patient offset percentile fractions less than threshold data.

        Parameters
        ----------
        offset_type : str, optional
            The motion metric percentiles to evaluate against thresholds. Default is
            "target_offset_vector".
        threshold_step : float, optional
            Step width to use to define thresholds against which to evaluate offset
            percentiles, in mm. Default is 1.
        figsize : tuple, optional
            Figure size. Default is (12, 4), i.e., 12 by 4 inches.

        Returns
        -------
        mpl.Figure
            Figure showing fraction of patient offset percentiles less than a
            threshold.

        Notes
        -----
        This figure is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        df = self.session_fraction_less_than_threshold(offset_type, threshold_step)
        return self._plot_fraction_less_than_threshold(
            df, "patients", offset_type, figsize
        )

    def plot_session_fraction_less_than_threshold(
        self,
        offset_type: str = "target_offset_vector",
        threshold_step: float = 1,
        figsize=(12, 4),
    ) -> mpl.Figure:
        """Plots session offset percentile fractions less than threshold data.

        Parameters
        ----------
        offset_type : str, optional
            _description_, by default "target_offset_vector"
        threshold_step : float, optional
            _description_, Default is 1 mm.
        figsize : tuple, optional
            Figure size. Default is (12, 4), i.e., 12 by 4 inches.

        Returns
        -------
        mpl.Figure
            Figure showing fraction of session offset percentiles less than a
            threshold.

        Notes
        -----
        This figure is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        df = self.session_fraction_less_than_threshold(offset_type, threshold_step)
        return self._plot_fraction_less_than_threshold(
            df, "sessions", offset_type, figsize
        )

    def session_fraction_less_than_threshold(
        self, offset_type: str = "target_offset_vector", threshold_step: float = 1
    ) -> pl.DataFrame:
        """Calculates the fraction of session motion statistics less than thresholds,
        intended to provide insight on margin selection.

        Parameters
        ----------
        offset_type : str, optional
            The motion metric percentiles to evaluate against thresholds. Default is
            "target_offset_vector".
        threshold_step : float, optional
            Step width to use to define thresholds against which to evaluate offset
            percentiles, in mm. Default is 1.

        Returns
        -------
        pl.DataFrame
            DataFrame containing defining fraction of session offset percentiles less
            than calculated thresholds.

        Notes
        -----
        This calculation is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        metrics = self.metrics.filter(pl.col("session_index").is_not_null())
        return self._fraction_less_than_threshold(metrics, offset_type, threshold_step)

    def to_csv(self, path: str | os.PathLike) -> None:
        """Writes motion data to CSV file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for CSV file to be written.

        Returns
        -------
        None
        """
        self._df.write_csv(path)

    def to_npz(self, path: str | os.PathLike, compress: bool = True) -> None:
        """Writes motion data to numpy npz file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for npz file to be written.
        compress : bool, optional
            Flag to indicate whether to compress npz. Default is True.
        """
        data_dict = {}
        if "timestamp" in self._df:
            data_dict["timestamp"] = self._df["timestamp"].to_numpy()
        else:
            data_dict["delta_time"] = self._df["delta_time"].to_numpy()
        for col in [
            "target_offset_x",
            "target_offset_y",
            "target_offset_z",
            "rigid_body",
            "potential_diff",
            "measured_diff",
            "potential_diff_tolerance",
            "measured_diff_tolerance",
        ]:
            if col in self._df:
                data_dict[col] = self._df[col].to_numpy()
        method = np.savez_compressed if compress else np.savez
        method(path, **data_dict)

    def to_parquet(self, path: str | os.PathLike) -> None:
        """Writes motion data to parquet file.

        Parameters
        ----------
        path : str | os.PathLike
            Path for parquet file to be written.
        """
        self._df.write_parquet(path)

    # endregion

    # region Private methods

    # endregion

    # region Static methods

    @staticmethod
    def _calculate_additional_offsets(df: pl.DataFrame) -> pl.DataFrame:
        """Calculates vector and delta offset values for dataframe.

        Parameters
        ----------
        df : pl.DataFrame
            DataFrame containing target offsets in x, y and z dimensions.

        Returns
        -------
        pl.DataFrame
            DataFrame containing vector and delta offset values.
        """
        result_df = df.with_columns(
            ((pl.col("target_offset_x") ** 2) + (pl.col("target_offset_z") ** 2))
            .sqrt()
            .alias("target_offset_xz"),
            (
                (pl.col("target_offset_x") ** 2)
                + (pl.col("target_offset_y") ** 2)
                + (pl.col("target_offset_z") ** 2)
            )
            .sqrt()
            .alias("target_offset_vector"),
            (pl.col("target_offset_x") - pl.col("target_offset_x").first()).alias(
                "delta_target_offset_x"
            ),
            (pl.col("target_offset_y") - pl.col("target_offset_y").first()).alias(
                "delta_target_offset_y"
            ),
            (pl.col("target_offset_z") - pl.col("target_offset_z").first()).alias(
                "delta_target_offset_z"
            ),
        )
        result_df = result_df.with_columns(
            (
                (pl.col("delta_target_offset_x") ** 2)
                + (pl.col("delta_target_offset_y") ** 2)
                + (pl.col("delta_target_offset_z") ** 2)
            )
            .sqrt()
            .alias("delta_target_offset_vector")
        )
        return result_df

    @staticmethod
    def _calculate_additional_times(df: pl.DataFrame) -> pl.DataFrame:
        """Calculates delta time and datetimes for dataframe.

        Parameters
        ----------
        df : pl.DataFrame
            DataFrame containing timestamps.

        Returns
        -------
        pl.DataFrame
            DataFrame containing delta time and datetimes.
        """
        result_df = df
        result_df = result_df.with_columns(
            ((pl.col("timestamp") - pl.col("timestamp").first()) / 1000).alias(
                "delta_time"
            )
        )
        result_df = result_df.with_columns(
            pl.from_epoch("timestamp", time_unit="ms").alias("datetime")
        )
        return result_df

    @staticmethod
    def _calculate_metrics(df: pl.DataFrame) -> pl.DataFrame:
        df_results = []
        # Calculate duration (or skip, if motion contains multiple deliveries)
        if (
            len(df.select("delta_time").to_numpy())
            == np.count_nonzero(df.select("delta_time").to_numpy()) + 1
        ):
            df_results.append(
                df.select(
                    [
                        (
                            pl.col("delta_time").last() - pl.col("delta_time").first()
                        ).alias("duration"),
                        pl.col("delta_time").count().alias("num_data_points"),
                    ]
                )
            )
        else:
            df_results.append(
                df.select(
                    [
                        pl.lit(None).alias("duration"),
                        pl.col("delta_time").count().alias("num_data_points"),
                    ]
                )
            )
        # Calculate IEC-X, IEC-Y, IEC-Z metrics
        for col in ["target_offset_x", "target_offset_y", "target_offset_z"]:
            df_results.append(
                df.select(
                    [
                        pl.col(col).mean().alias("mean_" + col),
                        pl.col(col).std().alias("standard_deviation_" + col),
                        pl.col(col).median().alias("median_" + col),
                        (pl.col(col).quantile(0.75) - pl.col(col).quantile(0.25)).alias(
                            "interquartile_range_" + col
                        ),
                        (pl.col(col).quantile(0.9) - pl.col(col).quantile(0.1)).alias(
                            "interdecile_range_" + col
                        ),
                        (pl.col(col).quantile(0.95) - pl.col(col).quantile(0.05)).alias(
                            "5pc_trimmed_range_" + col
                        ),
                        (pl.col(col).max() - pl.col(col).min()).alias("range_" + col),
                        pl.col("delta_" + col).mean().alias("mean_delta_" + col),
                        pl.col("delta_" + col).median().alias("median_delta_" + col),
                    ]
                )
            )
        # Calculate vector displacement and tracking metrics
        for col in [
            series
            for series in [
                "target_offset_vector",
                "delta_target_offset_vector",
                "potential_diff",
                "rigid_body",
                "measured_diff",
            ]
            if series in df.columns
        ]:
            df_results.append(
                df.select(
                    [
                        pl.col(col).mean().alias("mean_" + col),
                        pl.col(col).std().alias("standard_deviation_" + col),
                        pl.col(col).median().alias("median_" + col),
                        pl.col(col).quantile(0.8).alias("80th_percentile_" + col),
                        pl.col(col).quantile(0.9).alias("90th_percentile_" + col),
                        pl.col(col).quantile(0.95).alias("95th_percentile_" + col),
                        pl.col(col).max().alias("maximum_" + col),
                    ]
                )
            )
        # Calculate adaptation metrics
        df_results.append(
            df.select(
                [
                    (
                        (pl.col("target_offset_y").abs() > 12.5).sum()
                        / pl.col("delta_time").count()
                    ).alias("fraction_target_offset_y_exceeding_12.5"),
                    (
                        (pl.col("target_offset_y").abs() > 20).sum()
                        / pl.col("delta_time").count()
                    ).alias("fraction_target_offset_y_exceeding_20"),
                    (
                        (
                            (
                                (3.125 / pl.col("target_offset_xz"))
                                .arccos()
                                .fill_nan(0)
                                * 2
                                / np.pi
                            ).mean()
                        ).alias("fraction_target_offset_xz_exceeding_3.125")
                    ),
                    (
                        (
                            (
                                (4 / pl.col("target_offset_xz")).arccos().fill_nan(0)
                                * 2
                                / np.pi
                            ).mean()
                        ).alias("fraction_target_offset_xz_exceeding_4")
                    ),
                ]
            )
        )
        # Return combined calculated metrics
        return pl.concat(df_results, how="horizontal")

    @staticmethod
    def _fraction_less_than_threshold(
        metrics: pl.DataFrame,
        offset_type: str = "target_offset_vector",
        threshold_step: float = 1,
    ) -> pl.DataFrame:
        """Calculates the fraction of motion statistics less than thresholds,
        intended to provide insight on margin selection.

        Parameters
        ----------
        metrics : pl.DataFrame
            The motion metrics containing offset percentiles r80, r90, r95 and r100.
        offset_type : str, optional
            The motion metric percentiles to evaluate against thresholds. Default is
            "target_offset_vector".
        threshold_step : float, optional
            Step width to use to define thresholds against which to evaluate offset
            percentiles, in mm. Default is 1.

        Returns
        -------
        pl.DataFrame
            DataFrame containing defining fraction of offset percentiles less than
            calculated thresholds.

        Notes
        -----
        This calculation is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        r80 = metrics.select(pl.col(f"80th_percentile_{offset_type}")).to_numpy()
        r90 = metrics.select(pl.col(f"90th_percentile_{offset_type}")).to_numpy()
        r95 = metrics.select(pl.col(f"95th_percentile_{offset_type}")).to_numpy()
        r100 = metrics.select(pl.col(f"maximum_{offset_type}")).to_numpy()
        result_dict = {}
        thresholds = np.arange(
            0, int(np.max(np.ceil(r100))) + threshold_step, threshold_step
        )
        result_dict["thresholds"] = thresholds
        result_dict["fraction_r80"] = [
            np.count_nonzero(r80 < threshold) / len(r80) for threshold in thresholds
        ]
        result_dict["fraction_r90"] = [
            np.count_nonzero(r90 < threshold) / len(r90) for threshold in thresholds
        ]
        result_dict["fraction_r95"] = [
            np.count_nonzero(r95 < threshold) / len(r95) for threshold in thresholds
        ]
        result_dict["fraction_r100"] = [
            np.count_nonzero(r100 < threshold) / len(r100) for threshold in thresholds
        ]
        return pl.DataFrame(result_dict)

    @staticmethod
    def _plot_fraction_less_than_threshold(
        fraction_less_than_threshold_df: pl.DataFrame,
        grouping: str,
        offset_type: str,
        figsize=(12, 4),
    ) -> mpl.Figure:
        """Plots offset percentile fractions less than threshold data.

        Parameters
        ----------
        fraction_less_than_threshold_df : pl.DataFrame
            DataFrame containing offset percentile fraction less than threshold data.
        grouping : str, optional
            Description of grouping used for dataframe. E.g., "patients" or "session".
        offset_type : str, optional
            The motion metric percentiles to evaluate against thresholds.
        figsize : tuple, optional
            Figure size. Default is (12, 4), i.e., 12 by 4 inches.

        Returns
        -------
        mpl.Figure
            Figure showing fraction of offset percentiles less than a threshold.

        Notes
        -----
        This figure is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        offset_type_label = (
            (" ".join([offset_type.split("_")[-1]] + offset_type.split("_")[0:-1]))
            .replace("x ", "IEC-X ")
            .replace("y ", "IEC-Y ")
            .replace("z ", "IEC-Z ")
            .replace("delta ", "Relative ")
            .replace("vector ", "Vector ")
            .replace(" Vector", " vector")
        )
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.axhline(0.95, color="grey", linestyle="--")
        ax.plot(
            fraction_less_than_threshold_df["thresholds"],
            fraction_less_than_threshold_df["fraction_r80"],
            label="for 80% tracking time",
        )
        ax.plot(
            fraction_less_than_threshold_df["thresholds"],
            fraction_less_than_threshold_df["fraction_r90"],
            label="for 90% tracking time",
        )
        ax.plot(
            fraction_less_than_threshold_df["thresholds"],
            fraction_less_than_threshold_df["fraction_r95"],
            label="for 95% tracking time",
        )
        ax.plot(
            fraction_less_than_threshold_df["thresholds"],
            fraction_less_than_threshold_df["fraction_r100"],
            label="for 100% tracking time",
        )
        ax.set_xlim(
            fraction_less_than_threshold_df["thresholds"][0],
            fraction_less_than_threshold_df["thresholds"][-1],
        )
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))
        ax.set_xlabel(f"{offset_type_label} threshold (mm)")
        y_label = offset_type_label.replace("Relative", "relative").replace(
            "Vector", "vector"
        )
        ax.set_ylabel(f"Percentage of {grouping} with \n{y_label} less than threshold")
        ax.legend()
        return fig

    # endregion
