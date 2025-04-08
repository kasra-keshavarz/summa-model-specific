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
    Union,
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
        # Crude check for special occasion of exactextractr (MAF) or QGIS
        # provided inputs
        if isinstance(stats, pd.DataFrame):
            self.data = stats

        # Otherwise, check if `stats` is a Python dictionary
        if not isinstance(stats, pd.DataFrame) and not isinstance(stats, dict):
            raise TypeError("`stats` must be of dtype dict.")

        # Otherwise, all good and build the object and use pandas
        # DataFrame object for ease of use and flexibility
        self.data = pd.DataFrame.from_dict(stats)

        # `None` default value for self.threshold
        self.threshold = None

        # if `frac` is included, sort the fractions
        if any(self.data.columns.str.startswith('frac')):
            pass # do sorting

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
        return [Stats.__regex_extract_frac(i) for i in self['frac']]

    # static methods
    @staticmethod
    def __regex_extract_frac(
        string,
    ) -> int | str | float:
        """Extract values after frac_ values

        Parameters
        ----------
        string : str
           String starting with `frac_` based on the standards of this object.
        """
        return re.search(r"frac_([^_]+)", string).group(1)

    @staticmethod
    def __convert_list_dtype(
        ls,
    ) -> List[int | float | str]:
        """Converts str elements of a list to int or float dtypes, if
        feasible"""
        # Just for reassurance
        ls = list(ls)

        # Convert to `int` or `float` dtype if feasible
        converted = []

        # Start conversion
        for item in ls:
            try:
                # Try converting to int first
                converted.append(int(item))
            except ValueError:
                try:
                    # If int fails, try float
                    converted.append(float(item))
                except ValueError:
                    # If both fail, keep original
                    converted.append(item)
        return converted

    @staticmethod
    def mask_classes(
        row: pd.Series,
        mask_series: pd.Series
    ) -> float:
        """
        Return masked values based on `mask_series` where corresponding row values are non-zero.

        For each row in a DataFrame, masks the mask_series based on zero values in the row,
        then returns the sum of the remaining mask_series values.

        Parameters
        ----------
        row : pd.Series
            A single row from the DataFrame to be processed. The index should align with
            mask_series index.
        mask_series : pd.Series
            Series containing values to be summed after masking. Must share index with
            the DataFrame columns.

        Returns
        -------
        pd.Series
            The mask_series values where the corresponding row values were non-zero.

        Notes
        -----
        - The function is designed to be used with pandas.DataFrame.apply()
        - Zero values in the row will mask (exclude) corresponding values in mask_series
        - Only indices present in both the row and mask_series will be considered

        Examples
        --------
        >>> df = pd.DataFrame({'A': [1, 0], 'B': [0, 1]})
        >>> mask_series = pd.Series([10, 20], index=['A', 'B'])
        >>> df.apply(mask_classes, axis=1, mask_series=mask_series)
        0    10.0  # Only A=1 kept (10)
        1    20.0  # Only B=1 kept (20)
        dtype: float64
        """
        # Create a mask where row values are not zero
        mask = (row != 0)

        # Align the mask with mask_series (since they share the same index)
        aligned_mask = mask[mask_series.index]
 
        # Apply the mask to mask_series and sum the remaining values
        return mask_series[aligned_mask]

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

    # core functions
    def map_fracs(
        self: Self,
        mapping: Dict[Any, Any],
        inplace: bool = True,
        autoupdate_stats: bool = False,
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
        autoupdate_stats : bool [defaults to `False`]
            Autoupdate other statistics, if applicable, based on new `frac`
            values.

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
        # Assuring there are no `None` values - they must be in the `removing`
        # dictionary instead
        renaming = {k: v for k, v in renaming.items() if v is not None}
        # remove any `None` keys
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
        sorted_cols = [Stats.__regex_extract_frac(c) for c in df.columns]
        # Convert dtype from str to int or float if possible
        converted_cols = Stats.__convert_list_dtype(sorted_cols)
        # Sort the values
        converted_cols = [_prefix+str(c) for c in sorted(converted_cols)]
        # Sort the columns
        df = df[converted_cols]

        # After all the operations, if there are duplicate columns, raise an
        # Exception
        if len(df.columns) != len(set(df.columns)):
            raise ValueError("Duplicates cannot exist with `frac` stats.")

        # If some `classes` were removed, re-normalize fractions
        if len(removing.values()) >= 1:
            df = df.div(df.sum(axis=1), axis=0)

        # If certain rows turn to have np.nan values due to adjustments above,
        # turn them into zeros
        df.fillna(0, inplace=True)

        # Change ['frac'] values to `df`
        if inplace == True:
            self['frac'] = df
            return df
        else:
            return df
            
    def _update_stats(
        self: Self,
    ) -> None:
        """Adjust relevant statistics after `frac` changes"""
        # Due to removal of certain classes, we need to adjust the statistics
        # while also warning user of possible misalignment of stats

        # Empty set to keep track of what stats have been updated
        _updated_stats = {'frac'} 

        # Extract _classes to calculate new `mean` values - just as a
        # precautionary measure not to mess things up; otherwise, self.classes
        # is available
        _classes = Stats.__convert_list_dtype([Stats.__regex_extract_frac(c)
                                                  for c in self['frac'].columns]) 
        # Convert `_classes` to pandas.Series to ease multiplications
        _classes = pd.Series(_classes, index=self['frac'].columns)

        # Checking stats that needs to be updated
        if 'minority' in self.stats:
            _updated_stats.add('minority')
            _minority = self['frac'].replace(0, np.nan).idxmin(axis=1)

        if 'majority' in self.stats:
            _updated_stats.add('majority')
            _majority = self['frac'].replace(0, np.nan).idxmax(axis=1)

        if 'min' in self.stats:
            _updated_stats.add('min')
            # Selecting the first non-zero class for the element;
            # If element does not cover any classes (frac for all is 0), then
            # report `np.nan` and a warning
            min_lambda = lambda x: x.replace(0, np.nan).dropna().index[0] \
                if not x.replace(0, np.nan).dropna().empty else np.nan
            _min = self['frac'].apply(min_lambda, axis=1)
            # Send warning to user
            wrn_msg = ("Some elements do not cover any classes, therefore, "
                "corresponding `min` is set to `np.nan`.")
            warnings.warn(wrn_msg)

        if 'max' in self.stats:
            _updated_stats.add('max')
            # Selecting the last non-zero class for the element;
            # If element does not cover any classes (frac for all is 0), then
            # report `np.nan` and a warning
            max_lambda = lambda x: x.replace(0, np.nan).dropna().index[-1] \
                if not x.replace(0, np.nan).dropna().empty else np.nan
            _max = self['frac'].apply(max_lambda, axis=1)
            # Send warning to user
            wrn_msg = ("Some elements do not cover any classes, therefore, "
                "corresponding `max` is set to `np.nan`.")
            warnings.warn(wrn_msg)

        if 'mean' in self.stats:
            _updated_stats.add('mean')
            # Calculate new mean values
            _means = self['frac'].multiply(_classes).sum(axis=1)

        if 'variety' in self.stats:
            _updated_stats.add('variety')
            # Calculate a rough count based on `frac` values
            _variety = self['frac'].apply(Stats.mask_classes, axis=1, mask_series=_classes).apply(lambda s: s.count(), axis=1)
            variety_lambda = lambda x: x.replace(0, np.nan).dropna().index[-1]
            _variety = self['frac'].apply(variety_lambda, axis=1)
            # Send warning to user
            wrn_msg = ("Some elements do not cover any classes, therefore, "
                "corresponding `max` is set to `np.nan`.")
            warnings.warn(wrn_msg)

        # Check stats that are dependent on `count`
        if 'count' in self.stats:
            _updated_stats.add('count')
            # Calculate new `count`

            if 'q' in self.stats:
                _updated_stats.add('q')
                pass

            if 'median' in self.stats:
                _updated_stats.add('median')
                pass

            if 'stdev' in self.stats:
                _updated_stats.add('stdev')
                pass

            if 'coefficient_of_variation' in self.stats:
                _updated_stats.add('coefficient_of_variation')
                pass

            if 'variance' in self.stats:
                _updated_stats.add('variance')
                pass

        _left_out = self.stats - _updated_stats
        _left_out -= {'frac'}
        if _left_out:
            # Warn users about stats that were not autoupdated
            wrn_msg = f"The following stats were not updated: {_left_out}"
            warnings.warn(wrn_msg)

        # FIXME: new stat updates to be added later.

        return

    def set_frac_threshold(
        self: Self,
        threshold: float | int,
        inplace: bool = True,
        autoupdate_stats: bool = False,
    ) -> pd.DataFrame:
        """Setting the `frac_threshold` value which removes those `classes`
        having fractions less than `threshold` on average over all elements.

        threshold : float
            Threshold defining minimum fractions for 
        """
        # Check the dtype
        if not isinstance(threshold, int | float):
            raise TypeError("`threshold` ")
        # Checking whether the relevant statistics exist in the first place
        if 'frac' in self.stats:
            # Setting the value
            if inplace:
                self.threshold = threshold

            # Find average fraction over all elements
            _average = self['frac'].mean(axis=0).copy()

            # Find `frac` `class`es where they are less than `threshold`
            _threshold_classes = _average.index[_average >= threshold]

            # Resemble a mapping dictionary
            _threshold_classes = [Stats.__regex_extract_frac(c) 
                                      for c in _threshold_classes]
            _psuedo_mapping = {k: None for k in _threshold_classes}

            # cal map_fracs method
            df = self.map_fracs(_psuedo_mapping, inplace=inplace)

            return df

        else:
            raise ValueError("`frac_threshold` cannot be set without `frac` stats.")

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
        stats: Dict[str, Dict[str, Any]],
        layer: np.ndarray = None,
        geom: gpd.GeoDataFrame = None,
        engine: str = 'gdal',
    ) -> None:
        """Main constructor for GeoSpatial layers

        Parameters
        ----------
        layer : array-like
            An array-like object representing raster layer of interest
        geom : |GeoDataFrame|
            A GeoDataFrame representing elements for which 
        """
        super().__init__(
            stats=stats
        )
        # type of engine must be `str`
        if not isinstance(engine, str):
            raise ValueError("`engine` must have a dtype of `str`.")

        # assign attributes if provided
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

