"""RadixactDatasetCohort module.

This module provides functionality for processing cohorts of data from Radixact
treatments. Uses lazy loading to avoid too many datasets being opened concurrently.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import logging
import os
from functools import cached_property
from typing import Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from pyradlib.radixact.dataset import RadixactDataset
from pyradlib.radixact.motion import RadixactSynchronyMotion
from pyradlib.radixact.record import RadixactRecord

logger = logging.getLogger(__name__)


class RadixactDatasetCohort:
    # region Constructors

    def __init__(self, df: pl.DataFrame) -> RadixactDatasetCohort:
        """Initialises an object corresponding to a cohort of Radixact patients.

        Parameters
        ----------
        df : pl.DataFrame
            DataFrame containing the patients and corresponding paths comprising the
            cohort. Columns should contain "patient" for patient identifier, and "path"
            for location of patient dataset.

        Returns
        -------
        RadixactDatasetCohort
            The encapsulated cohort dataset object.
        """
        self._df = df
        logger.info(f"Cohort initialised with {len(self._df)} patients")

    @classmethod
    def from_path(cls, path: str | os.PathLike) -> RadixactDatasetCohort:
        """_summary_

        Parameters
        ----------
        path : str | os.PathLike
            Path to directory containing subdirectories corresponding to patient
            datasets.

        Returns
        -------
        RadixactDatasetCohort
            The encapsulated cohort dataset object.
        """
        return cls.from_path_list(
            [
                os.path.join(path, dir)
                for dir in sorted(os.listdir(path))
                if os.path.isdir(os.path.join(path, dir))
            ]
        )

    @classmethod
    def from_path_list(
        cls, paths: list[str] | list[os.PathLike]
    ) -> RadixactDatasetCohort:
        """Initialises an object corresponding to a cohort of Radixact patients.

        Parameters
        ----------
        paths : list[str] | list[os.PathLike]
            Paths to folders containing patient datasets.

        Returns
        -------
        RadixactDatasetCohort
            The encapsulated cohort dataset object.
        """
        result_dict = {}
        result_dict["patient"] = [os.path.basename(path.rstrip("\\")) for path in paths]
        result_dict["path"] = paths
        return cls(pl.DataFrame(result_dict))

    # endregion

    # region Magic methods

    def __len__(self) -> int:
        """Calculates the number of patients in the cohort as length.

        Returns
        -------
        int
            Length of cohort data, or number of patients.
        """
        return len(self._df)

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactDatasetCohort(num_datasets={len(self._df)}"

    # endregion

    # region Properties

    @cached_property
    def motion(self) -> RadixactSynchronyMotion:
        """Returns concatenated motion data.

        Returns
        -------
        RadixactSynchronyMotion
            Concated Synchrony motion data from the cohort.
        """
        if len(self._df) == 0:
            return None
        else:
            return RadixactSynchronyMotion.from_patient_motions(self.motions)

    @cached_property
    def motion_metrics(self) -> pl.DataFrame:
        """Returns concatenated motion metrics for patient cohort.

        Returns
        -------
        pl.DataFrame
            Concatenated motion metrics for patient cohort.
        """
        if len(self._df) == 0:
            return None
        else:
            return self.motion.metrics

    @cached_property
    def motions(self) -> list[RadixactSynchronyMotion]:
        """Returns list of motions for patient cohort.

        Returns
        -------
        list[RadixactSynchronyMotion]
            List of encapsulated motions, one for each patient.

        Note
        ----
        This method is slow. It requires parsing of as many XML files as there are
        delivery sessions in the patient cohort.
        """
        motions = []
        for dir in self._df["path"]:
            ds = RadixactDataset.from_path(dir)
            motions.append(ds.motion)
        return motions

    @cached_property
    def plan_summary(self) -> pl.DataFrame:
        """Returns summary of plans for patient cohort.

        Returns
        -------
        pl.DataFrame
            Summary of plans for patient cohort.
        """
        if len(self._df) == 0:
            return None
        else:
            plan_summaries = []
            for patient_index, dir in enumerate(self._df["path"]):
                ds = RadixactDataset.from_path(dir)
                plan_summaries.append(
                    ds.plan_summary.with_columns(patient_index=pl.lit(patient_index))
                )
            return pl.concat(plan_summaries).select(
                [pl.col("patient_index"), pl.all().exclude("patient_index")]
            )

    @cached_property
    def plan_settings_summary(self) -> pl.DataFrame:
        """Returns concatenated plan settings for patient cohort.

        Returns
        -------
        pl.DataFrame
            Concatenated plan settings for patient cohort.
        """
        if len(self._df) == 0:
            return None
        else:
            plan_settings_summaries = []
            for patient_index, dir in enumerate(self._df["path"]):
                ds = RadixactDataset.from_path(dir)
                plan_settings_summaries.append(
                    ds.plan_settings_summary.with_columns(
                        patient_index=pl.lit(patient_index)
                    )
                )
            return pl.concat(plan_settings_summaries).select(
                [pl.col("patient_index"), pl.all().exclude("patient_index")]
            )

    @cached_property
    def records_summary(self) -> pl.DataFrame:
        """Produce summary of treatment record parameters.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment record parameters.
        """
        if len(self._df) == 0:
            return None
        else:
            records_summaries = []
            for patient_index, dir in enumerate(self._df["path"]):
                ds = RadixactDataset.from_path(dir)
                records_summaries.append(
                    ds.records_summary.with_columns(patient_index=pl.lit(patient_index))
                )
            return pl.concat(records_summaries, how="vertical_relaxed").select(
                [pl.col("patient_index"), pl.all().exclude("patient_index")]
            )

    @cached_property
    def telemetry_metrics(self) -> pl.DataFrame:
        """Returns concatenated telemetry metrics for patient cohort.

        Returns
        -------
        pl.DataFrame
            Concatenated telemetric metrics.
        """
        if len(self._df) == 0:
            return None
        else:
            telemetry_metrics = []
            for patient_index, dir in enumerate(self._df["path"]):
                ds = RadixactDataset.from_path(dir)
                if len(ds.telemetry_sinograms) > 0:
                    telemetry_metrics.append(
                        ds.telemetry_metrics.with_columns(
                            patient_index=pl.lit(patient_index)
                        )
                    )
            return pl.concat(telemetry_metrics).select(
                [pl.col("patient_index"), pl.all().exclude("patient_index")]
            )

    # endregion

    # region Public methods

    def imaging_angles(self) -> np.ndarray:
        """Extracts imaging angles from plan settings summaries for each patient.

        Returns
        -------
        np.ndarray
            Array containing valid imaging angles.
        """
        result = (
            self.plan_settings_summary.select(
                [
                    "imaging_angle_1",
                    "imaging_angle_2",
                    "imaging_angle_3",
                    "imaging_angle_4",
                    "imaging_angle_5",
                    "imaging_angle_6",
                ]
            )
            .to_numpy()
            .flatten()
        )
        return result[~np.isnan(result)]

    def plot_correction_histogram(
        self,
        parameters: list[str] | None = None,
        binwidth: float = 0.25,
        col_wrap: int = 4,
        sharex: bool = False,
        sharey: bool = False,
    ) -> mpl.Figure:
        """Generate a histogram plot of treatment couch and gantry corrections after
        daily imaging and registration.

        Parameters
        ----------
        parameters : list[str], optional
            The couch adjustments to include in the figure. Possible options include:
            "correction_x", "correction_y", "correction_z"
                Couch adjustment in the IEC-X, IEC-Y and/or IEC-Z dimensions.
            "correction_roll"
                Gantry angle adjustment.
            Default is None, in which case correction_x, correction_y, correction_z,
            and correction_roll are used.
        binwidth: float, optional.
            Width of bin, in mm or degrees. Default is 0.25 (mm or degrees).
        col_wrap : int, optional
            Number of column facets within a row. Default is 4.
        sharex : bool, optional
            Flag for whether facets should share x axes. Default is False.
        sharey : bool, optional
            Flag for whether facets should share y axes. Default is True.

        Returns
        -------
        mpl.Figure
            Histogram of target offset values in each dimension.
        """
        return RadixactRecord.plot_correction_histogram(
            self.records_summary, parameters, binwidth, col_wrap, sharex, sharey
        )

    def plot_imaging_angles(self, width: int = 6) -> mpl.Figure:
        """_summary_

        Parameters
        ----------
        width : int, optional
            Width of plot bars in degrees. Should be an even positive factor of 360,
            or 1. Default is 6.

        Returns
        -------
        mpl.Figure
            Polar bar plot of radiographic imaging angles.
        """
        imaging_angles = self.imaging_angles()
        theta = np.linspace(0.0, 360, int(360 / width), endpoint=False)
        if width == 1:
            radii, _ = np.histogram(imaging_angles, 360, (0, 360))
        if width > 1:
            offset_angles = ((imaging_angles + int(width / 2)) % 360) - int(width / 2)
            radii, _ = np.histogram(
                offset_angles, int(360 / width), (-int(width / 2), 360 - int(width / 2))
            )
        bar_width = width * np.pi / 180
        fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "polar"})
        ax.bar(theta * np.pi / 180, radii, width=bar_width, bottom=0)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        return fig

    def plot_motion_boxplot_sns(
        self,
        parameters: list[str] | None = None,
        aspect: int = 1,
        col_wrap: int = 4,
        sharey: bool = False,
    ) -> mpl.Figure:
        """Plot boxplot of motion and tracking parameters using Seaborn.

        Parameters
        ----------
        parameters : list[str], optional
            The Synchrony paramters to include in the figure. Possible options include:
            "target_offset_x", "target_offset_y", "target_offset_z"
                Target offset in IEC-X, IEC-Y and/or IEC-Z dimensions.
            "target_offset_vector"
                Target offset in 3D space.
            "delta_target_offset_x", "delta_target_offset_y", "delta_target_offset_z"
                Target offset relative to starting position in IEC-X, IEC-Y and/or
                IEC-Z dimensions.
            "delta_target_offset_vector"
                Target offset in 3D space relative to starting position.
            "potential_diff"
                Potential difference or uncertainty in target position calculcation.
            "measured_diff"
                Measured delta or difference between observed and calculated offset.
            "rigid_body"
                Rigid body value describing deformation of fiducial distribution.
            Default is None, in which case target_offset_x, target_offset_y,
            target_offset_z, and target_offset_vector are used.
        aspect : int, optional
            Aspect ratio of each facet. Higher aspect ratios may require change to col_wrap. Default is 1.
        col_wrap : int, optional
            Number of column facets within a row. Default is 4.
        sharey : bool, optional
            Flag for whether facets should share y axes. Default is False.

        Returns
        -------
        mpl.Figure
            Boxplot of values in each dimension.
        """
        return self.motion.plot_motion_boxplot_sns(parameters, aspect, col_wrap, sharey)

    def plot_motion_histogram(
        self,
        mode: Literal["absolute", "delta", "relative", "starting"] = "absolute",
        fig_size: tuple[float, float] = (12, 4),
        offset_lim: tuple[float, float] | None = None,
        offset_bin: float = 0.25,
        vector_lim: tuple[float, float] | None = None,
        vector_bin: float = 0.25,
        title: str | None = None,
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
        mlp.Figure
            Histogram of target offset values in each dimension.
        """
        return self.motion.plot_motion_histogram(
            mode, fig_size, offset_lim, offset_bin, vector_lim, vector_bin, title
        )

    def plot_motion_histogram_sns(
        self,
        parameters: list[str] | None = None,
        binwidth: float = 0.25,
        col_wrap: int = 4,
        sharex: bool = False,
        sharey: bool = True,
    ) -> mpl.Figure:
        """Generate a histogram plot of distribution of offset values in each dimension,
        using the Seaborn module.

        Parameters
        ----------
        parameters : list[str], optional
            The Synchrony paramters to include in the figure. Possible options include:
            "target_offset_x", "target_offset_y", "target_offset_z"
                Target offset in IEC-X, IEC-Y and/or IEC-Z dimensions.
            "target_offset_vector"
                Target offset in 3D space.
            "delta_target_offset_x", "delta_target_offset_y", "delta_target_offset_z"
                Target offset relative to starting position in IEC-X, IEC-Y and/or
                IEC-Z dimensions.
            "delta_target_offset_vector"
                Target offset in 3D space relative to starting position.
            "potential_diff"
                Potential difference or uncertainty in target position calculcation.
            "measured_diff"
                Measured delta or difference between observed and calculated offset.
            "rigid_body"
                Rigid body value describing deformation of fiducial distribution.
            Default is None, in which case target_offset_x, target_offset_y,
            target_offset_z, and target_offset_vector are used.
        binwidth: float, optional.
            Width of bin, in mm. Default is 0.25 (mm).
        col_wrap : int, optional
            Number of column facets within a row. Default is 4.
        sharex : bool, optional
            Flag for whether facets should share x axes. Default is False.
        sharey : bool, optional
            Flag for whether facets should share y axes. Default is True.

        Returns
        -------
        mpl.Figure
            Histogram of target offset values in each dimension.
        """
        return self.motion.plot_motion_histogram_sns(
            parameters, binwidth, col_wrap, sharex, sharey
        )

    def plot_patient_motions_fraction_less_than_threshold(
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
        return self.motion.plot_patient_fraction_less_than_threshold(
            offset_type, threshold_step, figsize
        )

    def plot_session_motions_fraction_less_than_threshold(
        self,
        offset_type: str = "target_offset_vector",
        threshold_step: float = 1,
        figsize=(12, 4),
    ) -> mpl.Figure:
        """Plots the fraction of session motion statistics less than thresholds,
        intended to provide insight on margin selection.

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
            Figure showing fraction of session motion statistics less than thresholds.

        Notes
        -----
        This calculation is inspired by Figure 5(b) of Li et al. (2008), available at
        DOI:10.1016/j.ijrobp.2007.10.049.
        """
        return self.motion.plot_session_fraction_less_than_threshold(
            offset_type, threshold_step, figsize
        )

    def plot_percentile_vector_target_offset(
        self, fig_size: tuple[float, float] = (12, 4)
    ) -> mpl.Figure:
        """Calculates R95, R90, and R80 percentile vector target offsets for individual
        patients within the cohort, and presents as a line plot.

        Returns
        -------
        mpl.Figure
            Line plot for percentile vector target offset.

        Notes
        -----
        This plot resembles that of Figure 5(b) from Li et al.'s 2007 paper, Dosimetric
        consequences of intrafraction prostate motion.
        """
        patient_metrics = self.motion_metrics.filter(
            pl.col("session_index").is_null(), pl.col("patient_index").is_not_null()
        )
        fig, ax = plt.subplots(1, 1, figsize=fig_size)
        r80 = patient_metrics.select(
            pl.col("80th_percentile_target_offset_vector")
        ).to_numpy()
        r90 = patient_metrics.select(
            pl.col("90th_percentile_target_offset_vector")
        ).to_numpy()
        r95 = patient_metrics.select(
            pl.col("95th_percentile_target_offset_vector")
        ).to_numpy()
        pt = patient_metrics.select(pl.col("patient_index")).to_numpy() + 1
        ax.plot(pt, r80, label="$r_{80}$")
        ax.plot(pt, r90, label="$r_{90}$")
        ax.plot(pt, r95, label="$r_{95}$")
        ax.set_xlim(pt[0], pt[-1])
        ax.set_xlabel("Patient number")
        ax.set_ylabel("Vector target offset (mm)")
        ax.legend()
        return fig

    # endregion
