"""Geospatial related workflows
"""
# 3rd party libraries
import pandas as pd
import numpy as np

# built-in libraries
import os
from typing import (
    Self,
    Mapping,
    Any,
)


class Stats(object):
    """Stats used in the workflows and relevant methods defined

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(
        self: Self,
        stats: Dict[str, Dict[str, Any]],
    ) -> Self:
        """Main constructor for `Stats`

        Parameters
        ----------
        stats : dict
            Dictionary reporting one or more statistics for elements

        Returns
        -------
        |DataFrame|
            A pandas.DataFrame containing statistics of interest

        Notes
        -----
        - Acceptable keys for `stats` are:
          * **frac**
          * **mean**
          * **min**
          * **max**
          * **mean**
          * **majority**
          * **minority**
          * **median**
          * **quantile**
          * **variety**
          * **variance**
          * **stdev**
          * **coefficient_of_variation**
          * **frac**
          * **coords**
          * **count**
          * **sum**
        """
        # special occasion of MAF provided inputs
        if isinstance(stats, pd.DataFrame):
            self.stats = stats
            return 

        # otherwise, check if `stats` is a Python dictionary
        if not isinstance(stats, dict):
            raise TypeError("`stats` must be of dtype dict.")

        # otherwise, all good and build the object
        self.stats = pd.DataFrame.from_dict()

        return

    @property
    def 


    @classmethod
    def from_maf(
        gistool_csv: str | os.PathLike,
    ) -> Self:
        obj = pd.read_csv(gistool-csv,
            index_col=0,
            header=0)
        # special occasion for qgis and MAF ecists in __init__
        # this needs better coding, but I don't have enough time
        return cls(obj)

    @classmethod
    def from_qgis(
        qgis_csv: str | os.PathLike,
    ) -> Self:
        obj = pd.read_csv(qgis_csv,
            index_col=0,
            header=0)
        # special occasion for qgis and MAF exists in __init__
        return cls(obj)


class GeoLayer(object):
    """GeoLayer defining geospatial data
    """

    def __init__(
        self: Self,
        layer: np.ndarray = None,
        stats: Stats = None,
        engine: str = 'gdal',
    ) -> None:
        """Main constructor for GeoSpatial layers
        """
        # at least one input option should be provided
        if layer is None and stats is None:
            raise ValueError("Either `layer` or `stats` should be provided.")

        # type of engine must be `str`
        if not isinstance(engine, str):
            raise ValueError("`engine` must have a dtype of `str`.")

        # assign attributes if provided
        if stats:
            self.stats = Stats
        if layer:
            if engine.lower() == 'gdal'
            self.layer = layer

        return

    def __repr__(
        self
    ) -> str:
        """Official representation of the object
        """
        # Show the stats provided
        stats_provided = list(self.stats.keys())
        # Show the layer value
        layer_provided = 'True' if self.layer else 'False'

        # reprs
        stats_repr = f"Stats: {stats_provided}"
        layer_repr = f"Layer: {layer_provided}"

        # return the string representation of the layer
        return stats_repr + '\n' + layer_repr

    def plot(
        self
    ) -> None:
        """Plots the `self.layer`"""
        try:
           import matplotlib.pyplot as plt
        except ImportError:
           raise ImportError("To plot GeoLayer object, you need to install matplotlib>=3.1")

        if self.engine.lower() in ('gdal'):
            plt.imshow(self.layer)
            try:
                plt.show()
            else:
                raise ValueError("Cannot plot in non-interactive environment")

        return
