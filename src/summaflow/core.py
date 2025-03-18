"""
This package automates SUMMA model setup through a flexible workflow that can
be set up using a JSON configuration file, a Command Line Interface (CLI),
directly inside a Python script/environment, or other commonly used
configuration interfaces

The workflow is adopted from the efforts of Wouter Knoben, Darri
Eythorsson, and Mohamed Ismaiel Ahmed at the University of Calgary,
Canada.
"""
import pandas as pd
import numpy as np
import geopandas as gpd
import xarray as xr
import pint

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
    ddb_global_attrs_default,
    ddb_local_attrs_default,
    forcing_local_attrs_default,
    forcing_global_attrs_default,
    default_attrs,
)

class SUMMAWorkflow(object):
    """
    Main workflow class of SUMMA

    Attributes
    ----------

    Parameters
    ----------

    """
    # main constructor
    def __init__(
            self,
            riv: str | os.PathLike | gpd.GeoDataFrame,
            cat: str | os.PathLike | gpd.GeoDataFrame,
            landcover: str | os.PathLike,
            forcing_files: Sequence[str | os.PathLike],
            forcing_vars: Sequence[str],
            main_id: str,
            ds_main_id: str,
            landcover_classes: Dict[str, str] = None,
            forcing_units: Dict[str, str] = None,
            forcing_to_units: Dict[str, str] = None,
            outlet_value: int = -9999,
            riv_cols: Dict[str, Union[str, int]] = None,
            cat_cols: Dict[str, Union[str, int]] = None,
            local_attrs: Dict[str, Dict] = None,
            global_attrs: Dict[str, Dict] = None,
            gru_dim: str = 'gru',
            hru_dim: str = 'subbasin',
            ) -> None:

        """Main constructor of SUMMAWorkflow
        """
