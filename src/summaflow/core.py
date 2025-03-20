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
    forcing_data : a sequence of :obj:`str` or |PathLike|
        A sequence to NetCDF files describing forcing data for the
        SUMMA application of interest.
    forcing_vars : :obj:`dict` of :obj:`str`
        The keys are the variable names included in `forcing_files`
        that should be included in the final forcing dataset(s). The
        corresponding values are the accepted names SUMMA needs.
    forcing_units : :obj:`dict` of :obj:`str` keys and values
        The keys are the **accepted** variable names for SUMMA, and
        values are the original units provided in `forcing_files`.
        The units follow the standards of Pint's `default_en.txt`
        manual. See the relevant reference.
    forcing_to_units : :obj:`dict` of :obj:`str` keys and values
        The keys are the **accepted** variable names for SUMMA, and
        values are the **accepted**  units for SUMMA. The units
        follow the standards of Pint's `default_en.txt` manual. See
        the relevant reference.
    topology_data : :obj:`dict` of :obj:`str` sequences
        Topology data generally consists of auxillary data needed to
        analyze `riv`, `cat`, and `hru` objects. The keys are the
        strings of each mentioned objects, and values include data
        labels in each object needed to be included in the final
        model setup.
    topology_units : :obj:`dict`
        Original units for all data labels included in `topology_data`
        elements
    topology_to_units : :obj:`dict`
        **Acceptable** units of all data labels included in
        `topology_data` elements SUMMA requires
    geospatial_data : :obj:
        


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
        forcing_data: Sequence[str | os.PathLike],
        forcing_attrs: Dict[str] = None,
        topology_data: Dict[str, str | int] = None,
        topology_attrs: Dict[str, str] = None,
        geospatial_data: Dict[str, Dict] = None,
        geospatial_attrs: Dict[str, Dict] = None,
        dims: Dict[str, str] = None,
    ) -> None:

        """
        Main constructor of SUMMAWorkflow
        """
