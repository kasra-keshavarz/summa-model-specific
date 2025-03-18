"""
This package automates SUMMA model setup through a flexible workflow that can
be set up using a JSON configuration file, a Command Line Interface (CLI),
directly inside a Python script/environment, or other commonly used
configuration interfaces
"""
import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import pint
import jinja2 as jij

from typing import (
    Dict,
    Sequence,
    Union,
)

# built-in libraries
import re
import json
import sys
import glob
import os
import shutil
import warnings

from ._default_dicts import (
    forcing_local_attrs_default,
    forcing_global_attrs_default,
    default_attrs,
)

class SUMMAWorkflow(object):
    """
    Main workflow class of SUMMA

    Attributes
    ----------
    riv : :obj:`str` or |PathLike| or |GeoDataFrame|
        River graph object. Either the path to the object should be
        given or a geopandas.GeoDataFrame containing information
        provided in `topology`. This aids in understanding the
        topology of the river system and hierarchical order of its
        connectivity.
    cat : :obj:`str` or |PathLike| or |GeoDataFrame|
        Basin polygon(s) object. Either the path to the object
        should be given or a geopandas.GeoDataFrame containing
        information provided in `topology`. This object usually
        serves as GRU information for SUMMA applications.
    hru : :obj:`str` or |PathLike| or |GeoDataFrame|
        HRU polygon(s) object. Either the path to the object
        should be given or a geopandas.GeoDataFrame containing
        information provided in `topology`. This object serves
        as HRU information for SUMMA applications. If it equals to
        `None`, then each HRU will be equal to GRU elements.
    forcing_files : a sequence of :obj:`str` or |PathLike|
        A sequence to NetCDF files describing forcing data for the
        SUMMA application of interest.
    forcing_vars : :obj:`dict` of :obj:`str`
        The keys are the variable names included in `forcing_files`
        that should be included in the final forcing dataset(s). The
        corresponding values are the accepted names SUMMA needs.
    forcing_units : :obj:`dict` of :obj:`str`
        The keys are the **accepted** variable names for SUMMA, and
        values are the original units provided in `forcing_files`.
        The units follow the standards of Pint's `default_en.txt`
        manual. See the relevant reference.
    forcing_to_units : :obj:`dict` of :obj:`str`
        The keys are the **accepted** variable names for SUMMA, and
        values are the **accepted**  units for SUMMA. The units
        follow the standards of Pint's `default_en.txt` manual. See
        the relevant reference.

    Properties
    ----------


    References
    ----------
    .. [1] `pint/default_en.txt <https://github.com/hgrecco/pint/blob/master/pint/default_en.txt>`_
           The default unit definitions used by the pint package.

    See Also
    --------
    pint : The pint package for handling units in Python.

    """
    # main constructor
    def __init__(
        self,
        riv: str | os.PathLike | gpd.GeoDataFrame,
        cat: str | os.PathLike | gpd.GeoDataFrame,
        hru: str | os.PathLike | gdp.GeoDataFrame,
        forcing_files: Sequence[str | os.PathLike],
        forcing_vars: Sequence[str],
        forcing_units: Dict[str, str] = None,
        forcing_to_units: Dict[str, str] = None,
        topology: Dict[str, str | int] = None,
        topology_attrs: Dict[str, str] = None,
        topology_units: Dict[str, str] = None,
        topology_to_units: Dict[str, str] = None,
        local_attrs: Dict[str, Dict] = None,
        global_attrs: Dict[str, Dict] = None,
        dims: Dict[str, str] = None,
    ) -> None:

        """
        Main constructor of SUMMAWorkflow

        Note
        ----
        `hru` needs to provide mapping information to `cat`. Please note
        that `hru` can also be equal to `cat`.
        
        outlet_value: int = -9999,
        main_id: str,
        ds_main_id: str,
        riv_cols: Dict[str, Union[str, int]] = None,
        cat_cols: Dict[str, Union[str, int]] = None,
        gru_dim: str = 'gru',
        hru_dim: str = 'hru',
        """
