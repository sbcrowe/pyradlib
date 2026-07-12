# -*- coding: utf-8 -*-
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

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from pyradlib.radixact.dataset import RadixactDataset
from pyradlib.radixact.motion import RadixactSynchronyMotion

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
        logger.debug(f"Cohort contains {len(self._df)} patients")

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
                for dir in os.path.sorted(os.listdir(path))
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
        return f"RadixactDatasetCohort(num_datasets={str(len(self._df))}"

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
        mlp.Figure
            Histogram of target offset values in each dimension.
        """
        return self.motion.plot_motion_histogram(
            mode, fig_size, offset_lim, offset_bin, vector_lim, vector_bin, title
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
