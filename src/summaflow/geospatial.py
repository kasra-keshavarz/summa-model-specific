"""Geospatial related workflows"""
# built-in libraries
import os
import re
import warnings
from collections import (
    OrderedDict,
)
from collections.abc import (
    Sequence,
)
from typing import (
    Any,
    Dict,
    List,
    Optional,
)

try:
    from typing import Self  # Python 3.11+
except ImportError:
    from typing_extensions import Self  # For Python <3.11

# 3rd party libraries
import geopandas as gpd
import numpy as np
import pandas as pd
import pint
import pint_pandas

# import internal functions
from . import utils

# constant values and functions
## pandas.DataFrame index slicer
idx = pd.IndexSlice
## NoneType is only defined on Python@3.10+
NoneType = type(None)
## List of acceptable statistics for the objects defined
ACCEPTABLE_STATS = {
    'min',
    'max',
    'mean',
    'majority',
    'minority',
    'median',
    'quantile',
    'variety',
    'variance',
    'stdev',
    'coefficient_of_variation',
    'frac',
    'coords',
    'count',
    'sum',
}
## List of acceptable statistics accepting units
STATS_WITH_UNITS = {
    'min',
    'max',
    'mean',
    'majority',
    'minority',
    'median',
    'quantile',
    'sum',
    'coords',
}
## List of acceptable statistics with squared units
STATS_WITH_SQUARED_UNITS = {
    'variance',
}
## List of acceptable statistics that are dimensionless
STATS_DIMENSIONLESS = {
    'stdev',
    'frac',
    'count',
    'coefficient_of_variation',
}
STATS_WITHOUT_UNITS = STATS_DIMENSIONLESS


class Stats:
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

    # virtual properties
    @property
    def stats_info(
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
            # Since `frac_%` can be weird, remove all from `self.stats_info`,
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
        string: str,
    ) -> int | str | float:
        """Extract values after frac_ values

        Parameters
        ----------
        string : str
           String starting with `frac_` based on the standards of this object.

        Notes
        -----
        - Name mangling is implement to avoid confusion with subclass objects.
        """
        return re.search(r"frac_([^_]+)", string).group(1)

    @staticmethod
    def __convert_list_dtype(
        ls: List[Any],
    ) -> List[int | float | str]:
        """Converts str elements of a list to int or float dtypes, if
        feasible

        Notes
        -----
        - Name mangling is implement to avoid confusion with subclass objects.
        """
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
        - Name mangling is implement to avoid confusion with subclass objects.

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

        # If more than one argument is entered, throw an IndexError
        if isinstance(item, tuple) and len(item) > 1:
            raise IndexError("Too many indices provided.")

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

        # If `item` not in `self.stats_info` raise IndexError
        for stat in items:
            if stat not in self.stats_info:
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
        view = self.data.loc[:, cols].squeeze()

        if view.shape == (1, 1):
            view = view.iloc[:, 0]
        else:
            view = view.squeeze()

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
        # if `stat` is in `self.stats_info` remove them, if not available proceed
        # with inserting
        if stat in self.stats_info:
            # define relevant columns to be removed
            cols = self[stat].columns
            # drop the relevant columns
            self.data.drop(columns=cols, inplace=True, errors='raise')

        # add the relevant columns of stats to the `.data`
        self.data = pd.concat([self.data, df], axis=1)


    def __repr__(
        self: Self,
    ) -> pd.DataFrame:
        """Official representation of the object"""
        number_of_stats = len(self.stats_info)
        return f"{number_of_stats} statistic(s) provided"

    # core functions
    def map_fracs(
        self: Self,
        mapping: Dict[Any, Any],
        inplace: Optional[bool] = False,
        autoupdate_stats: Optional[bool] = False,
    ) -> Self:
        """Implement fraction renaming, summation, and renormalization
        of values, based on `mapping`

        Parameters
        ----------
        mapping : dict
            Mapping values from the class fractions already available to a
            target class system.
        inplace : bool [defaults to `False`]
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
        if 'frac' not in self.stats_info:
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

        # Update object's stats
        if inplace == True:
            # note that other stats, in case of `autoupdate_stats=True` are
            # taken care of with `self._update_stats(...)` function
            self['frac'] = df

        # Autoupdate other stats if asked
        if autoupdate_stats:
            df = pd.concat(
                        [df, self._update_stats(fracs=df, inplace=inplace)],
                        axis=1)

        # just return a copy of the caluclated adjustments
        return df

    def _update_stats(
        self: Self,
        fracs: pd.DataFrame,
        inplace: Optional[bool] = False,
    ) -> None:
        """Adjust relevant statistics after `frac` changes.
        
        This method updates various statistics (like minority, majority, min,
        max, etc.) based on the current `frac` values. It handles cases where
        classes may have been removed and warns about potential misalignments
        in statistics.

        Parameters
        ----------
        fracs : |DataFrame|
            New `frac` database, based on which the rest of `self.stats_info`
            elements are updated.
        inplace : bool [defaults to `False`]
            Change the stats in place.

        Notes
        -----
        - Updates statistics that are dependent on `frac` values
        - Handles edge cases where elements don't cover any classes (setting
          to np.nan)
        - Warns about statistics that couldn't be automatically updated
        - Uses lambda functions for efficient computation of certain
          statistics
 
        Warnings
        --------
        - Issues warnings when:
            * Elements don't cover any classes (min/max/variety become np.nan)
            * Certain statistics couldn't be updated automatically
            * Potential misalignment between statistics

        Examples
        --------
        >>> stats_obj._update_stats()
        ... # Updates internal statistics based on current frac values

        The following statistics may be updated if present in self.stats_info:
        - minority: class with smallest non-zero fraction
        - majority: class with largest non-zero fraction
        - min: first non-zero class
        - max: last non-zero class
        - mean: weighted mean of classes
        - variety: count of non-zero classes
        """
        # Due to removal of certain classes, we need to adjust the statistics
        # while also warning user of possible misalignment of stats

        # Set to keep track of what stats have been updated
        # Since this function may be called after changing `frac` values, the
        # set contains it to begin with.
        _updated_stats = {'frac'}

        # Extract _classes to calculate new `mean` values - just as a
        # precautionary measure not to mess things up; otherwise, self.classes
        # is available
        _classes = Stats.__convert_list_dtype([Stats.__regex_extract_frac(c)
                                                  for c in self['frac'].columns])
        # Convert `_classes` to pandas.Series to ease multiplications
        _classes = pd.Series(_classes, index=self['frac'].columns)

        # Lambdas
        lambdas = {
            'minority': lambda x: Stats.__regex_extract_frac(x.replace(0, np.nan).dropna().idxmin()) \
                if not x.replace(0, np.nan).dropna().empty else np.nan,
            'majority': lambda x: Stats.__regex_extract_frac(x.replace(0, np.nan).dropna().idxmax()) \
                if not x.replace(0, np.nan).dropna().empty else np.nan,
            'min': lambda x: Stats.__regex_extract_frac(x.replace(0, np.nan).dropna().index[0]) \
                if not x.replace(0, np.nan).dropna().empty else np.nan,
            'max': lambda x: Stats.__regex_extract_frac(x.replace(0, np.nan).dropna().index[-1]) \
                if not x.replace(0, np.nan).dropna().empty else np.nan,
            'variety': lambda x: x.replace(0, np.nan).dropna().count(),
            'mean': None, # for compeleteness
            'count': None, # FIXME: to be completed
            'q': None, # FIXME: to be completed
            'stdev': None, # FIXME: to be completed
            'variance': None, # FIXME: to be completed
            'coefficient_of_variation': None, # FIXME: to be completed
            'median': None, # FIXME: to be completed
        }

        # If `inplace` is True
        data = self.data.copy()

        # Checking stats that needs to be updated
        if 'minority' in self.stats_info:
            _updated_stats.add('minority')
            data['minority'] = fracs.apply(lambdas['minority'], axis=1)

        if 'majority' in self.stats_info:
            _updated_stats.add('majority')
            data['majority'] = fracs.apply(lambdas['majority'], axis=1)

        if 'min' in self.stats_info:
            _updated_stats.add('min')

            # Selecting the first non-zero class for the element;
            # If element does not cover any classes (frac for all is 0), then
            # report `np.nan` and a warning
            data['min'] = fracs.apply(lambdas['min'], axis=1)

            # Send warning to user
            wrn_msg = ("Some elements do not cover any classes, therefore, "
                "corresponding `min` is set to `np.nan`.")
            warnings.warn(wrn_msg)

        if 'max' in self.stats_info:
            _updated_stats.add('max')

            # Selecting the last non-zero class for the element;
            # If element does not cover any classes (frac for all is 0), then
            # report `np.nan` and a warning
            data['max'] = fracs.apply(lambdas['max'], axis=1)

            # Send warning to user
            wrn_msg = ("Some elements do not cover any classes, therefore, "
                "corresponding `max` is set to `np.nan`.")
            warnings.warn(wrn_msg)

        if 'mean' in self.stats_info:
            _updated_stats.add('mean')

            # Calculate new mean values
            data['mean'] = fracs.multiply(_classes).sum(axis=1)

        if 'variety' in self.stats_info:
            _updated_stats.add('variety')

            # Calculate `variety` based on adjusted `frac` values
            data['variety'] = fracs.apply(lambdas['variety'], axis=1)
            # _variety = self['frac'].apply(lambdas['variety'], axis=1) # FIXME: needs to be addressed

            # Send warning to user
            wrn_msg = ("Some elements do not cover any classes, therefore, "
                "corresponding `max` is set to `np.nan`.")
            warnings.warn(wrn_msg)

        # Check stats that are dependent on `count`
        if 'count' in self.stats_info:
            # _updated_stats.add('count')
            # Calculate new `count`

            if 'q' in self.stats_info:
                #_updated_stats.add('q')
                pass

            if 'median' in self.stats_info:
                #_updated_stats.add('median')
                pass

            if 'stdev' in self.stats_info:
                #_updated_stats.add('stdev')
                pass

            if 'coefficient_of_variation' in self.stats_info:
                #_updated_stats.add('coefficient_of_variation')
                pass

            if 'variance' in self.stats_info:
                #_updated_stats.add('variance')
                pass

        _left_out = self.stats_info - _updated_stats
        if _left_out:
            # Warn users about stats that were not autoupdated
            wrn_msg = f"The following stats were not updated: {_left_out}"
            warnings.warn(wrn_msg)

            # Add left over stats
            for stat in list(_left_out):
                data.update(self[stat])

        if inplace:
            self.data = data

        return data

    def set_frac_threshold(
        self: Self,
        threshold: float | int,
        inplace: Optional[bool] = False,
        autoupdate_stats: Optional[bool] = False,
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
        if 'frac' in self.stats_info:
            # Setting the value
            if inplace:
                # Just an informative variable for users; not used anywhere
                # else
                self.threshold = threshold

            # Find average fraction over all elements
            _average = self['frac'].mean(axis=0).copy()

            # Find `frac` `class`es where they are less than `threshold`
            _threshold_classes = _average.index[_average >= threshold]

            # Resemble a mapping dictionary
            _threshold_classes = [Stats.__regex_extract_frac(c)
                                      for c in _threshold_classes]
            _psuedo_mapping = {k: None for k in _threshold_classes}

            # Call map_fracs method
            df_fracs = self.map_fracs(_psuedo_mapping, inplace=inplace)

            # Autoupdate other stats if asked
            if autoupdate_stats:
                df_stats = self._update_stats(fracs=df_fracs, inplace=inplace)
                df = pd.concat([df_fracs, df_stats], axis=1)

            else:
                df = df_fracs

            return df

        else:
            raise ValueError("`frac_threshold` cannot be set without `frac` stats.")


class GeoLayer:
    """GeoLayer defining geospatial data as a subclass of Stats

    Attributes
    ----------

    Methods
    -------
    """
    def __init__(
        self: Self,
        stats: Stats,
        layer: Optional[np.ndarray] = None,
        geom: Optional[gpd.GeoDataFrame] = None,
        engine: Optional[str] = 'gdal',
        unit: Optional[str] = None,
    ) -> None:
        """Main constructor for GeoSpatial layers

        Parameters
        ----------
        layer : array-like
            An array-like object representing raster layer of interest
        geom : |GeoDataFrame|
            A GeoDataFrame representing elements for which stats are
            caluclated
        engine : str

        unit : str or None

        Notes
        -----
        - for dimensionless layers, use 'dimensionless' as unit, not None
        """
        # If `stats` is already of dtype `Stats` ignore building it;
        # This usually happens when one of the alternative constructors are
        # used.

        # stats attribute that refers to the superclass's Stats object
        self.stats = stats

        # Assign attributes if provided
        #if layer is not None:
        #    self.layer = layer
        self.layer = layer
        #if geom is not None:
        #    self.geom = geom
        self.geom = geom

        # GeoLayer engine for reading Geospatial files
        self._engine = engine

        # Pint unit registry
        self._ureg = pint.UnitRegistry()

        # if unit is `None`, warn users that this does not covert units
        if unit == None:
            warnings.warn("No unit conversion is taking place. See "
            "`.to_unit(...)` for such operation")

        # If nothing is given,
        # Assign first instance unit
        self._unit = self._ureg(unit).units


    # class methods
    @classmethod
    def from_maf(
        cls: Self,
        maf_stats: str | os.PathLike,
        maf_layer: Optional[str | os.PathLike] = None,
        maf_geolayer: Optional[str | os.PathLike] = None,
        *args,
        **kwargs,
    ) -> Self:
        """MAF-sepcific layer for internal users' convenience.

        Parameters
        ----------

        Returns
        -------
        GeoLayer
            A GeoLayer object based on MAF's outputs.
        """
        # Provide necessary values for superclass's `from_maf` constructor
        _stats_obj = Stats.from_maf(gistool_csv=maf_stats)

        # Populate `layer` and `geom` entities if provided (not `None`)
        return cls(
            stats=_stats_obj,
            layer=maf_layer,
            geom=maf_geolayer,
            *args,
            **kwargs)

    # special methods
    def __repr__(
        self: Self,
    ) -> str:
        """Official representation of the object
        """
        # __repr__ elements in an ordered dictionary for representation
        repr_dict = OrderedDict(
            stats_repr = f"Stats: {self.stats.stats_info}",
            layer_repr = f"Layer: {True if self.layer else False}",
            geom_repr = f"Geometry: {True if self.geom else False}",
            unit_repr = f"Geolayer Unit: {self._unit}",
        )

        # Main __repr__ string
        main_repr = "\n".join(repr_dict.values())

        # return the string representation of the layer
        return main_repr

    # static methods
    @staticmethod
    def assign_unit(
        series: pd.Series,
        unit: str | pint.Unit,
    ) -> pint_pandas.PintArray:
        """Apply a Pint unit to a |Series|"""
        # type check
        if not isinstance(unit, str | pint.Unit):
            raise TypeError("`unit` must be of type string or `pint.Unit`.")
        if not isinstance(series, pd.Series):
            raise TypeError("`col` must be of type pandas.Series")

        # If everything is OK
        series_unit = pint_pandas.PintArray(series, dtype=f"pint[{unit}]")

        return series_unit

    @staticmethod
    def convert_unit(
        series: pd.Series,
        target_unit: str,
    ) -> pd.Series:
        """Convert a |Series| unit to the `target_unit`."""
        # Routine checks
        if not isinstance(target_unit, str | pint.Unit):
            raise TypeError("`target_unit` must be of type string or `pint.Unit`.")
        if not isinstance(series, pd.Series):
            raise TypeError("`series` must be of type pandas.Series.")
        if not hasattr(series, 'pint'):
            raise ValueError("`series` must have a `pint` attribute (has a pint.Unit datatype).")

        # If everything is OK
        return series.pint.to(target_unit)

    @staticmethod
    def remove_unit(
        series: pd.Series,
    ) -> pd.Series:
        """Remove a pandas.Series unit value"""
        if not hasattr(series, 'pint'):
            raise ValueError("`series` must have a `pint` attribute (has a pint.Unit datatype).")

        return series.pint.magnitude

    # virtual properties (including getter and setter functions)
    @property
    def unit(
        self: Self,
    ) -> pint.Unit:
        """GeoLayer's Pint unit"""
        # report unit of the layer
        return self._unit

    @unit.setter
    def unit(
        self: Self,
        unit_value: str,
    ) -> None:
        """Setter function for .unit property"""
        # If unit_value is not string, throw a ValueError exception
        if not isinstance(unit_value, str):
            raise TypeError("`unit` must be a string value.")

        # If a unit is already assigned, you gotta make sure conversion takes
        # place

        # turn into a Pint.Unit object
        self._unit = self._ureg(unit_value).units


    @property
    def engine(
        self: Self,
    ) -> str:
        """Instance's `engine` property getter function"""

        return self._engine

    @engine.setter
    def engine(
        self: Self,
        _engine: str
    ) -> None:
        """Instance's `engine` property setter function"""
        if not isinstance(_engine, str):
            raise TypeError("`engine` must be a string value.")

        if _engine not in (""):
            raise ValueError("`engine` must be either `gdal` or `rasterio`.")

        self._engine = _engine


    # Object's methods
    def plot(
        self: Self,
        engine: str = 'matplotlib',
    ) -> None:
        """Plots the `self.layer`"""
        try:
           import matplotlib.pyplot as plt
        except ImportError:
           raise ImportError("To plot GeoLayer, install matplotlib>=3.1")

        if self.engine.lower() in ('gdal'):
            plt.imshow(self.layer)
            try:
                plt.show()
            except:
                raise ValueError("Cannot plot in non-interactive environment")

    def to_unit(
        self: Self,
        target_unit: str | NoneType,
        inplace: bool = True,
    ) -> None:
        """Tranfer GeoLayer's information to the new unit

        Parameters
        ----------
        target_unit : str
            To target unit to which the GeoLayer's information will be
            converted.
        inplace : bool [defaults to `True`]
            Changing values in place, if set to True, otherwise, return a copy
            of changes implemented.

        Notes
        -----
        - This modifies `unit` attribute in place.
        - Unit conversion logic is up to the end user.
        """
        if not isinstance(target_unit, str | NoneType | pint.Unit):
            raise TypeError("`unit` must be of type str, `None` "
            "(dimensionless), or `pint.Unit`.")

        # Change the stats values based on the new unit
        if self.stats:
            # Investigate the available stats and follow hard-coded rules for
            # each; This information needs to be hard-coded as computers do
            # not hold statistical knowledge; build unit mapping dictionary:
            unit_mapping = {}
            _temp_data = pd.DataFrame()

            # For each statistics, build the unit lambda function
            for stat in self.stats.stats_info:
                if stat in STATS_WITH_UNITS:
                    unit_mapping[stat] = lambda x: x
                if stat in STATS_WITHOUT_UNITS:
                    unit_mapping[stat] = lambda x: None
                if stat in STATS_WITH_SQUARED_UNITS:
                    unit_mapping[stat] = lambda x: x**2

                # Extract stats value for unit operations
                df = self.stats[stat].copy()

                # `stat.data` value
                if isinstance(df, pd.Series):
                    df = df.to_frame()
                elif isinstance(self.stats[stat], pd.DataFrame):
                    pass
                else:
                    raise TypeError("Stats type unacceptable.")

                # Assign units
                temp_df = df.apply(lambda col: GeoLayer.assign_unit(col, unit_mapping[col.name](self.unit)))
                # Implement unit conversion
                temp_df = temp_df.apply(lambda col: GeoLayer.convert_unit(col, unit_mapping[col.name](target_unit)))
                # Remove unit values
                temp_df = temp_df.apply(GeoLayer.remove_unit)

                # Add to _temp_data
                _temp_data = pd.concat([_temp_data, temp_df], axis=1)

        # Change the `stats.data` value
        if inplace:
            self.stats.data = _temp_data

        # Change the self.unit attribute
        self.unit = target_unit

        return _temp_data

