"""Geospatial related workflows"""
# 3rd party libraries
import pandas as pd
import numpy as np
import geopandas as gpd

# built-in libraries
import os
import re
import warnings
from typing import (
    Dict,
    Set,
    Self,
    Mapping,
    Any,
    List,
)
from collections.abc import (
    Sequence,
)
from collections import (
    Counter,
)

# import internal functions
from . import utils

# constant values and functions
idx = pd.IndexSlice


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
            Dictionary reporting one or more statistics for elements.

        Returns
        -------
        |DataFrame|
            A pandas.DataFrame containing statistics of interest.

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
          * **q%%**
          * **variety**
          * **variance**
          * **stdev**
          * **coefficient_of_variation**
          * **frac_%**
          * **coords**
          * **count**
          * **sum**

          The values are dictionaries with keys as the element for which the
          statistic is reported and value of the statistic of interest.

        - For **frac_%**, the suffix string explains the fraction of a
          quantity.
        - For **q%**, the suffix string denotes a quantile; e.g., 10th quantile
          is represented with `q10`.
        """
        # crude special occasion of exactextractr (MAF) or QGIS provided
        # inputs
        if isinstance(stats, pd.DataFrame):
            self.data = stats
            return 

        # otherwise, check if `stats` is a Python dictionary
        if not isinstance(stats, dict):
            raise TypeError("`stats` must be of dtype dict.")

        # otherwise, all good and build the object and use pandas
        # DataFrame object for ease of use and flexibility
        self.data = pd.DataFrame.from_dict(stats)

        return

    # virtual properties
    @property
    def stats(
        self,
    ) -> List[str]:
        """Items available in the object

        Returns
        -------
        Set : The available statistics in form of a set.
        """
        # Return a list of strings for statistics provided.
        # Order of items are immaterial, so using a Set.
        stats = set()

        # This can get illegible if written in a oneliner
        for s in self.data.columns:
            # Since `frac_%` can be weird, remove all from `self.stats`,
            # and only put `frac` there. Same applies to `q%`.
            if s.lower().startswith('frac'):
                item = 'frac'
            elif s.lower().startswith('q'):
                item = 'q'
            else:
                item = s
            # append to the stats Set
            stats.add(item)

        return stats

    @property
    def classes(
        self,
    ) -> List[int | str | float]:
        """Extract all known classes based on ['frac'] values"""
        return [Stats._regex_extract_frac(i) for i in self['frac']]

    # static methods
    @staticmethod
    def _regex_extract_frac(
        string,
    ) -> int | str | float:
        """Extract values after frac_ values

        Parameters
        ----------
        string : str
           String starting with `frac_` based on the standards of this object.
        """
        return re.search(r"frac_([^_]+)", string).group(1)

    # special methods
    def __getitem__(
        self,
        item: str | Sequence[str],
    ) -> pd.DataFrame:
        """Allowing statistics selection from the pandas.DataFrame object
        using []. This is needed, as `frac` and `q` can be difficult to deal
        with.

        Parameters
        ----------
        item : str or a Sequence of str
            List of statistics available in the object that needs to be
            viewed.

        Notes
        -----
        - FIXME: The custom entries for `frac` and `q` need improvement.
        """
        # If a sequence of str is given, make a tuple out of it
        if isinstance(item, Sequence):
            # make a set out of `item`
            items = set(item)

        # If only one element (string) is entered, make a tuple out of it
        if isinstance(item, str):
            # single element set
            items = {item,}

        # If not a string or Sequence, raise a TypeError
        if not isinstance(item, Sequence | str):
            raise TypeError('Item(s) must be a string or a Sequence of strings.')
 
        # Aliases for `frac` and `q`
        # For users with less knowledge of the whole thing
        aliases = {
            'frac': ('fractions', 'frac_', 'fracs', 'fraction', 'fracs_'),
            'q': ('quantiles', 'quantile'),
        }
        for stat, alias_options in aliases.items():
            if any(alias in items for alias in alias_options):
                # add proper stat name
                items.add(stat)
                # remove aliases if asked by the user
                items -= set(alias_options)

        # If `item` not in `self.stats` raise IndexError
        for stat in items:
            if stat not in self.stats:
                raise IndexError(f'`{stat}` cannot be indexed.')

        # Check the input item and see if it is available
        # Special case for `frac` and `q`
        cols = [] # empty list for convnience
        for stat in items:
           if stat in ("frac", "q"):
               cols += self.data.columns[self.data.columns.str.startswith(stat)].to_list()
           else:
               cols.append(stat)

        # returning a view not a copy
        view = self.data.loc[:, cols] 

        return view

    def __setitem__(
        self,
        stat: str,
        df: pd.DataFrame,
    ) -> None:
        """Adding statistics of interest to the object while updating other
        relevant statistics.

        Parameters
        ----------
        stat : |DataFrame|
           a new stats DataFrame to be added to the object

        Notes
        -----
        - 
        """
        # if `stat` is in `self.stats` remove them, if not available proceed
        # with inserting
        if stat in self.stats:
            # define relevant columns to be removed
            cols = self[stat].columns
            # drop the relevant columns
            self.data.drop(columns=cols, inplace=True, errors='raise')

        # add the relevant columns of stats to the `.data`
        self.data = pd.concat([self.data, df], axis=1)
 
        return

    def __repr__(
        self: Self,
    ) -> pd.DataFrame:
        """Official representation of the object"""
        number_of_stats = len(self.stats)
        return f"{number_of_stats} statistics provided"

    # custom constructors
    @classmethod
    def from_maf(
        cls: Self,
        gistool_csv: str | os.PathLike,
    ) -> Self:
        obj = pd.read_csv(gistool_csv, index_col=0, header=0)
        # special occasion for qgis and MAF exists in __init__
        # FIXME: this needs a bit more systematic approach
        return cls(obj)
    # alias for `from_maf` as it uses `exactextractr` under the hood
    from_extactextractr = from_maf

    @classmethod
    def from_qgis(
        qgis_csv: str | os.PathLike,
    ) -> Self:
        obj = pd.read_csv(qgis_csv, index_col=0, header=0)
        # special occasion for qgis and MAF exists in __init__
        return cls(obj)

    @classmethod
    def from_rasterstats(
        rasterstats_obj,
    ) -> Self:
        """Constructor for rasterstats Python package objects
        FIXME: This needs to be contributed by those who use the package
        """
        obj = rasterstats_obj # fix this
        return cls(obj)

    def map_fracs(
        self: Self,
        mapping: Dict[Any, Any],
        inplace: bool = True,
    ) -> Self:
        """Implement fraction renaming, summation, and renormalization
        of values, based on `mapping`

        Parameters
        ----------
        mapping : dict
            Mapping values from the class fractions already available to a
            target class system.
        inplace : bool [defaults to `True`]
            Doing the operations on the object. If not, set to `False`.

        Notes
        -----
        - If a value in `mapping` is used only once, simply the class will be
          renamed.
        - If a value is repeated multipled time for various keys, this means
          the classes represented in corresponding keys need to be summed up.
        - If a value of ``None`` is provided, the corresponding class key will
          be removed from the fraction classes, and all other values will be
          renormalized.
        - If a class is missing from the keys, that class will be assumed to
          have a ``None`` value.
        """
        # Constants
        _prefix = 'frac_'

        # Check the input data type
        if not isinstance(mapping, dict):
            raise TypeError("`mapping` must be of type dict.")

        # For mapping operations, `frac` stat is necessary
        if 'frac' not in self.stats:
            raise ValueError("`frac` values are necessary for mappping operations.")

        # Based on `mapping`, find the classes that should be renamed,
        # summed, and removed.
        # keys: `frac` classes to rename, values: final name
        renaming = utils.unique_dict_values(mapping)
        # keys: final name of the `classes` element, values: one or more
        # `classes` element
        summation = utils.nonunique_dict_values(mapping)
        # keys: `None`, values: all `classes` to be removed
        removing = {None: [k for k, v in mapping.items() if v is None]}

        # get a copy of .data to work on
        df = self[_prefix].copy()

        # First, do necessary `summation` for relevant fractions
        temp_df_summation = pd.DataFrame()
        for key, value in summation.items():
            # Create a pandas.DataFrame with _prefix+str(key) and summed value
            temp_df_summation[_prefix+str(key)] = df.loc[:, [_prefix + str(i) for i in value]].sum(axis=1).copy()

        # Second, delete `removing` fractions
        # Since the `.values()` of `removing` is always a list of `classes`
        # elements, with a single corresponding `.keys()` of `None`, we can
        # take the single element `.values()` and turn it into a list and call
        # its first ([0]) element.
        df.drop(
            columns=[_prefix+str(i) for i in list(removing.values())[0]],
            inplace=True,
            errors='raise',
        )

        # Third, build `renaming` fractions
        df.rename(columns={
                      _prefix+str(k): _prefix+str(v)
                      for k, v in renaming.items()},
                  inplace=True,
                  errors='raise',
        )

        # Finally, Remove the summation columns
        for key, value in summation.items():
            df.drop(
                columns=[_prefix+str(i) for i in value],
                inplace=True,
                errors='raise',
            )

        # Add the summed items back to the main `df`
        df = pd.concat([df, temp_df_summation], axis=1)
        # Sort the column names just for aesthetics
        df.sort_index(axis=1, inplace=True, ascending=False)

        # After all the operations, if there are duplicate columns, raise an
        # Exception
        if len(df.columns) != len(set(df.columns)):
            raise ValueError("Duplicates cannot exist with `frac` stats.")

        # If some `classes` were removed, re-normalize fractions
        if len(removing.values()) >= 1:
            df = df.div(df.sum(axis=1), axis=0)

        # Change ['frac'] values to `df`
        if inplace == True:
            self['frac'] = df
            return df
        else:
            return df

    @property
    def frac_threshold(
        self: Self,
    ) -> None:
        """If fraction threshold is defined or not"""
        return frac_threshold

    @frac_threshold.setter
    def _frac_threshold_setter(
        self: Self,
        threshold: float,
    ) -> None:
        """ """

        return

class GeoLayer(Stats):
    """GeoLayer defining geospatial data as a subclass of Stats

    Attributes
    ----------

    Methods
    -------
    """
    def __init__(
        self: Self,
        layer: np.ndarray = None,
        geom: gpd.GeoDataFrame = None,
        engine: str = 'gdal',
    ) -> None:
        """Main constructor for GeoSpatial layers

        Parameters
        ----------

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
            if engine.lower() == 'gdal':
                self.layer = layer

        return

    # special methods
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

    # plotting the `self.layer` object
    def plot(
        self: Self,
        engine: str = 'matplotlib',
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
            except:
                raise ValueError("Cannot plot in non-interactive environment")

        return
