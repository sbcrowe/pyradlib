"""Plan settings module.

This module provides functionality for processing of plan settings data from
Radixact treatments.
"""

# authorship information
__author__ = "Scott Crowe"
__email__ = "sb.crowe@gmail.com"
__credits__ = []
__license__ = "GPL3"

# import required code
import os
from functools import cached_property

import polars as pl


class RadixactPlanSettings:
    # region Constructors

    def __init__(self, df: pl.DataFrame) -> RadixactPlanSettings:
        """Initialise plan settings class.

        Parameters
        ----------
        df : pl.DataFrame
            The plan settings dataframe.

        Returns
        -------
        RadixactPlanSettings
            Plan settings encapsulated in helping wrapper.
        """
        self._df: pl.DataFrame = df

    @classmethod
    def from_plan_settings(cls, path: str | os.PathLike) -> RadixactPlanSettings:
        """Read plan settings from a plan.settings file.

        Parameters
        ----------
        path : path: str | os.PathLike
            Path to the plan.settings file.

        Returns
        -------
        RadixactPlanSettings
            Plan settings dataframe encapsulated in a RadixactPlanSettings wrapper.
        """
        data = {}
        with open(path, "r") as file:
            for line in file:
                if line.strip():
                    parameter, value = line.strip().split("=")
                    data[parameter.strip()] = (
                        value.strip().replace("false", "False").replace("true", "True")
                    )

        return cls(pl.DataFrame(data))

    # endregion

    # region Magic methods

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactPlanSettings(id={id(self)})"

    # endregion

    # region Properties

    @cached_property
    def summary(self) -> pl.DataFrame:
        """Produce summary of treatment plan settings.

        Returns
        -------
        pl.DataFrame
            DataFrame containing treatment plan settings.
        """
        result_dfs = []
        # Extract imaging protocols
        result_dfs.append(
            self._df.select(
                [
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(1)
                    .alias("kvct_protocol"),
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(2)
                    .alias("kvct_body_size"),
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(3)
                    .alias("kvct_field_of_view"),
                    pl.col("SharedCT-kvctKey")
                    .str.split(":")
                    .list.get(4)
                    .alias("kvct_mode"),
                    (
                        pl.col("scanSelectionMaxEdge").cast(pl.Float64)
                        - pl.col("scanSelectionMinEdge").cast(pl.Float64)
                    ).alias("scan_length"),
                ]
            )
        )
        # Extract Synchrony data, if available angles
        if "measuredDifferenceTolerance" in self._df:
            result_dfs.append(
                self._df.select(
                    [
                        pl.col("measuredDifferenceTolerance")
                        .cast(pl.Float64)
                        .alias("measured_diff_tolerance"),
                        pl.col("potentialDifferenceTolerance")
                        .cast(pl.Float64)
                        .alias("potential_diff_tolerance"),
                        pl.col("targetOffsetTolerance")
                        .cast(pl.Float64)
                        .alias("target_offset_tolerance"),
                        pl.col("rigidBodyError")
                        .cast(pl.Float64)
                        .alias("rigid_body_tolerance"),
                        pl.col("trackingRange")
                        .cast(pl.Float64)
                        .alias("tracking_range"),
                        (
                            pl.col("TrackingViolationDelayMicroSeconds").cast(
                                pl.Float64
                            )
                            / 1000000
                        ).alias("auto_pause_delay"),
                        pl.col("sensitivity"),
                    ]
                )
            )
            result_dfs.append(
                self._df.select(
                    pl.col(imaging_angle_name)
                    .cast(pl.Float64)
                    .alias(f"imaging_angle_{index + 1}")
                    if imaging_angle_name in self._df.columns
                    else pl.lit(None)
                    .cast(pl.Float64)
                    .alias(f"imaging_angle_{index + 1}")
                    for index, imaging_angle_name in enumerate(
                        [
                            "RadiographAngle_0-Angle",
                            "RadiographAngle_1-Angle",
                            "RadiographAngle_2-Angle",
                            "RadiographAngle_3-Angle",
                            "RadiographAngle_4-Angle",
                            "RadiographAngle_5-Angle",
                        ]
                    )
                )
            )
        return pl.concat(result_dfs, how="horizontal")

    # endregion

    def to_csv(self, path: str | os.PathLike):
        """Wrie plan settings to comma-separated value file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the CSV file.
        """
        self._df.to_csv(path)
