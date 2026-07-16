# -*- coding: utf-8 -*-
"""Timing analysis module.

This module provides functionality for processing of telemetry timing data from
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

import numpy as np
import polars as pl


class RadixactTiming:
    # region Constructors
    def __init__(self, df: pl.DataFrame, uid: str = None) -> RadixactTiming:
        """Initialises an object corresponding to telemetry timing.

        Parameters
        ----------
        df : DataFrame
            DataFrame containing telemetry timing data.
        uid : str, optional
            Unique identifier string, to allow association with telemetry fluence data.
            Default is None.

        Returns
        -------
        RadixactTiming
            DataFrame encapsulated with helper functions.
        """
        self._df = df
        self._uid = uid

    @classmethod
    def from_dat(cls, path: str | os.PathLike, uid: str = None) -> RadixactTiming:
        """Reads telemetry timing from a .dat file.

        Parameters
        ----------
        path : str | os.PathLike
            Path to the .dat file.
        uid : str, optional
            Unique identifier string, to allow association with telemetry fluence data.
            Default is None.

        Returns
        -------
        RadixactTiming
            The dat file is returned as a DataFrame encapsulated with helper functions.
        """
        with open(path) as timing_path_text_file:
            timing_text = timing_path_text_file.read()
        timing_text = timing_text.split("\n")
        # TODO Need to implement support for n > 1 data sets
        timestamps = []
        taus = []
        for n in range(int(timing_text[1])):
            # check for empty record
            if len(timing_text[3 + n * 9]) > 0:
                timing_timestamps = np.array(
                    [
                        x * 1000000 + y
                        for x, y in zip(
                            list(map(int, timing_text[3 + n * 9].split(","))),
                            list(map(int, timing_text[4 + n * 9].split(","))),
                        )
                        if x > 0
                    ]
                )
                timing_start_taus = np.array(
                    [
                        x
                        for x in list(map(float, timing_text[5 + n * 9].split(",")))
                        if x > 0
                    ]
                )
                timing_end_taus = np.array(
                    [
                        x
                        for x in list(map(float, timing_text[6 + n * 9].split(",")))
                        if x > 0
                    ]
                )
                timing_intervals = np.array(
                    [
                        x
                        for x in list(map(float, timing_text[7 + n * 9].split(",")))
                        if x > 0
                    ]
                )
                for timestamp, start_tau, end_tau, interval in zip(
                    timing_timestamps,
                    timing_start_taus,
                    timing_end_taus,
                    timing_intervals,
                ):
                    timestamps.append(timestamp)
                    timestamps.append(
                        timestamp + int((end_tau - start_tau) * interval * 1000000)
                    )
                    taus.append(start_tau)
                    taus.append(end_tau)
        data = pl.DataFrame({"tau": taus, "timestamp": timestamps})
        data = data.with_columns(
            pl.from_epoch("timestamp", time_unit="us").alias("datetime")
        )
        return cls(data, uid)

    # endregion

    # region Magic methods

    def __repr__(self) -> str:
        """Returns an unambigious string representation of the object.

        Returns
        -------
        str
            Representation of the object.
        """
        return f"RadixactTiming(uid={self._uid}, start={self._df['datetime'][0]})"

    # region Properties

    @cached_property
    def metrics(self) -> pl.DataFrame:
        """Calculates delivery timing metrics, describing treatment start, and total,
        beam, and interrupt durations.

        Returns
        -------
        pl.DataFrame
            The calculated timing metrics.
        """
        result = {}
        result["start_datetime"] = self._df.get_column("datetime").first()
        timestamps = self._df["timestamp"].to_numpy()
        result["total_duration"] = (timestamps[-1] - timestamps[0]) / 1000000
        paired_timestamps = timestamps.reshape(int(len(timestamps) / 2), 2)
        result["treatment_duration"] = (
            np.sum(paired_timestamps[:, 1] - paired_timestamps[:, 0]) / 1000000
        )
        result["interrupt_duration"] = (
            result["total_duration"] - result["treatment_duration"]
        )
        return pl.DataFrame(result)

    # endregion
