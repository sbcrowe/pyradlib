# -*- coding: utf-8 -*-
"""Motion analysis module.

This module provides functionality for processing of motion data from Synchrony
treatments.
"""

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
        x_displacements: npt.ArrayLike,
        y_displacements: npt.ArrayLike,
        z_displacements: npt.ArrayLike,
        rigid_body_values: npt.ArrayLike = None,
        potential_diffs: npt.ArrayLike = None,
        measured_deltas: npt.ArrayLike = None,
        uid: str = None,
    ) -> RadixactSynchronyMotion:
        """_summary_

        Parameters
        ----------
        times : npt.ArrayLike
            Relative timestamps of motion data points, in sec.
        x_displacements : npt.ArrayLike
            IEC-X or left-right displacements, in mm.
        y_displacements : npt.ArrayLike
            IEC-Y or superior-inferior displacements, in mm.
        z_displacements : npt.ArrayLike
            IEC-Z or anterior-posterior displacements, in mm.
        rigid_body_values : npt.ArrayLike, optional
            Rigid body fiducial deformation, in mm.
        potential_diffs : npt.ArrayLike, optional
            Potential difference of uncertainty in model, in mm.
        measured_deltas : npt.ArrayLike, optional
            Measured difference between model and images, in mm.
        uid : str, optional
            Unique identifier string, to allow association with other data. Default is
            None.

        Returns
        -------
        RadixactSynchronyMotion
            The Radixact Synchrony motion data encapsulated with helper functions.
        """
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
        df_dict["target_offset_x"] = x_displacements
        df_dict["target_offset_y"] = y_displacements
        df_dict["target_offset_z"] = z_displacements
        if rigid_body_values is not None:
            df_dict["rigid_body"] = rigid_body_values
        if potential_diffs is not None:
            df_dict["potential_diff"] = potential_diffs
        if measured_deltas is not None:
            df_dict["measured_diff"] = measured_deltas
        df = pl.DataFrame(df_dict)
        df = df.with_columns(
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
        df = df.with_columns(
            (
                (pl.col("target_offset_x") ** 2)
                + (pl.col("target_offset_y") ** 2)
                + (pl.col("target_offset_z") ** 2)
            )
            .sqrt()
            .alias("delta_target_offset_vector")
        )
        return cls(df, uid)

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
                    df = pl.DataFrame(data)
                    df = df.with_columns(
                        (
                            (pl.col("target_offset_x") ** 2)
                            + (pl.col("target_offset_z") ** 2)
                        )
                        .sqrt()
                        .alias("target_offset_lateral"),
                        (
                            (pl.col("target_offset_x") ** 2)
                            + (pl.col("target_offset_y") ** 2)
                            + (pl.col("target_offset_z") ** 2)
                        )
                        .sqrt()
                        .alias("target_offset_vector"),
                        (
                            pl.col("target_offset_x")
                            - pl.col("target_offset_x").first()
                        ).alias("delta_target_offset_x"),
                        (
                            pl.col("target_offset_y")
                            - pl.col("target_offset_y").first()
                        ).alias("delta_target_offset_y"),
                        (
                            pl.col("target_offset_z")
                            - pl.col("target_offset_z").first()
                        ).alias("delta_target_offset_z"),
                    )
                    df = df.with_columns(
                        (
                            (pl.col("target_offset_x") ** 2)
                            + (pl.col("target_offset_y") ** 2)
                            + (pl.col("target_offset_z") ** 2)
                        )
                        .sqrt()
                        .alias("delta_target_offset_vector")
                    )
                    dfs.append(df)
        # Create dataframe
        df = pl.concat(dfs)
        df = df.with_columns(
            ((pl.col("timestamp") - pl.col("timestamp").first()) / 1000).alias(
                "delta_time"
            )
        )
        df = df.with_columns(
            pl.from_epoch("timestamp", time_unit="ms").alias("datetime")
        )
        # Remove any non-existent data
        df = df.fill_nan(None)
        df_cleaned = df.select(
            [pl.col(col) for col in df.columns if not df[col].is_null().all()]
        )
        return cls(df_cleaned, uid)

    @classmethod
    def from_csv(
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
        df_results = []
        df_results.append(
            self._df.select(
                [
                    (pl.col("delta_time").last() - pl.col("delta_time").first()).alias(
                        "duration"
                    ),
                    pl.col("delta_time").count().alias("num_data_points"),
                ]
            )
        )
        # Calculate IEC-X, IEC-Y, IEC-Z metrics
        for col in ["target_offset_x", "target_offset_y", "target_offset_z"]:
            df_results.append(
                self._df.select(
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
            if series in self._df.columns
        ]:
            df_results.append(
                self._df.select(
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
            self._df.select(
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
                                (3.125 / pl.col("target_offset_lateral"))
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
                                (4 / pl.col("target_offset_lateral"))
                                .arccos()
                                .fill_nan(0)
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

    # endregion

    # region Public methods

    def to_csv(self, path) -> None:
        """Writes motion data to CSV file.

        Parameters
        ----------
        path : str
            Path for CSV file to be written.

        Returns
        -------
        None
        """
        self._df.write_csv(path)

    # endregion

    # region Private methods

    # endregion

    # region Static methods

    # endregion
