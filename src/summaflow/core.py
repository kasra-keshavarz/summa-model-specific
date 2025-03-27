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

# package imports
from ._default_dicts import (
    attributes_global_attrs_default,
    attributes_local_attrs_default,
    cold_state_local_attrs_default,
    cold_state_global_attrs_default,
    trial_params_local_attrs_default,
    trial_params_global_attrs_default,
    forcing_local_attrs_default,
    forcing_global_attrs_default,
    default_dims,
)
from .utils import (
    _init_empty_ds,
)


class SUMMAWorkflow(object):
    """
    Main workflow class of SUMMA

    Attributes
    ----------
    forcing : xr.DataSet
        Forcing dataset for the model


    Methods
    -------


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
        forcing_attrs: Dict[str, str],
        forcing_vars: Dict[str, str],
        forcing_units: Dict[str, str],
        forcing_to_units: Dict[str, str],
        topology_data: Dict[str, str | int],
        topology_vars: Dict,
        topology_units: Dict[str, str],
        topology_to_units: Dict[str, str],
        topology_attrs: Dict[str, str],
        geospatial_data: Dict[str, Dict] = None,
        dims: Dict[str, str] = default_dims,
    ) -> None:

        """
        Main constructor of SUMMAWorkflow

        Parameters
        ----------
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
        forcing_attrs : :obj:`dict`
            The keys are **acceptable** values of forcing variables
            SUMMA needs, and values are dictionaries of local attributes
            that need to be included in SUMMA's forcing file.
        topology_data : :obj:`dict` of :obj:`str` or |GeoDataFrame| Sequences
            Topology data generally consists of auxillary data needed to
            analyze `riv`, `cat`, and `hru` objects. The keys are the
            :obj:`str` of each mentioned objects, and values include data
            labels in each object needed to be included in the final
            model setup.
        topology_vars : :obj:`dict` of :obj:`str` keys and values
            Topology variables needed from the objects read through
            `topology_data`. Critical keys are: *segId*, *downSegId*,
            *slope*, *length*, *hruToSegId*, *tan_slope*. Any other
            additional labels will be added to SUMMA's Attribute NetCDF
            file for user's convenience.
        topology_units : :obj:`dict`
            Original units for all data labels included in `topology_data`
            elements
        topology_to_units : :obj:`dict`
            **Acceptable** units of all data labels included in
            `topology_data` elements SUMMA requires
        geospatial_data : :obj:`dict`
            Additional geospatial data used in parameterizing SUMMA.
            Acceptable keys are: *heightCanopyTop*
        dims : :obj:`dict`, optional
            Critical dimnesion names. Acceptable keys are *hru* and *gru*.

        Note
        ----
        The following lists the elements of `topology_data`
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
        `topology_attrs` needs the following keys:
           1. *gru_fid* that is a common label shared between `riv` and `cat`,
           2. *hru_fid* that maps `hru` values to `gru`,
        `forcing_attrs` needs the following keys:
           * mHeight

        Example
        -------

        """
        # routine error checks

        # assign necessary attributes
        # FIXME: This needs to turn into its own object, but for the
        #        sake of timing of this deliverable, we compromise and treat
        #        them as different variables starting with `forcing_`.
        self._forcing = forcing_data
        self._forcing_attrs = forcing_attrs
        self.forcing_vars = forcing_vars
        self.forcing_units = forcing_units
        self.forcing_to_units = forcing_to_units
        
        # topology data
        # `riv` is optional, as single-site configurations, do not have
        # any river systems involved
        try:
            self.riv = topology_data.get('riv')
        except:
            self.riv = None

        # `cat` will be added as a @property
        try:
            # `gru` can also be called `cat`
            self.gru = topology_data.get('cat')
        except:
            self.gru = topology_data.get('gru')

        # `hru` object
        self.hru = topology_data.get('hru')

        # `topology_attrs` object
        self.topology_attrs = topology_attrs

        # dimension names
        self.dims = dims

        # `init` object to add lazy-like behaviour
        self.init = False

        return


    @property
    def cat(
        self,
    ) -> gpd.GeoDataFrame:
        """aka `gru`"""
        return self.gru


    def init_attrs(
        self,
    ) -> None:
        """Initialize the necessary objects for the experient"""
        # variables to build an empty SUMMA-specific `attribute` object
        # since `self.topology_attrs['fid']` assumes to be found in both.
        # We use `self.gru` to build the `gru` variable values
        # the keys to `variables` is mandated by `utils._init_empty_ds`
        # function

        # FIXME: this is a terrible way of doing things, `gru` and `riv`
        #        should be taken care of by Hydrant, and `hru` should be
        #        taken care of using a separate descritization utility.
        #        The variable/dimension names should also be standards
        #        defined from upstream package(s).
        gru_fid = self.topology_attrs.get('gru_fid')
        hru_fid = self.topology_attrs.get('hru_fid')

        # defining `variables`
        variables = {
            'gru': self.gru[gru_fid],
            'hru': self.hru[hru_fid],
        }

        # create an empty xarray.Dataset to be populated with SUMMA-specific
        # attributes
        self.attrs = _init_empty_ds(
                         variables=variables,
                         dims=default_dims,
                     )

        # populating with necessary information obtained upon
        # instantiation; This class variable is a combination of
        # information of `forcing_*` and `topology_*` objects
        #
        # 1. measurement height value through self.forcing_attrs;
        #    the name for `mHeight` itself is taken from
        #    `forcing_vars['measurement_height']`. This value is defaulted to
        #    `mHeight` in the workflow, if it is not provided.
        # Notes:
        #    1. The `measurement_height` is the assumption of this workflow
        #    that is clearly instructed in the tutorial and API reference.
        #    2. `mHeight` is defined for each `hru` in SUMMA, therefore,
        #       hard-coded here.
        mHeight_name = self.forcing_vars.get('measurement_height')
        mHeight_values = [self._forcing_attrs.get('measurement_height')] * len(self.gru)
        self.attrs[mHeight_name] = xr.DataArray(mHeight_values, dims=self.dims['hru'])

        # 2. `slopeTypeIndex` which is "a legacy that is no longer used"
        #     reference: github.com/CH-Earth/CWARHM: step 5/SUMMA/1/1
        # Note:
        #     If needed, the hard-coded names & values here can be transformed
        #     to be read from the input objects.
        slope_type_index_values = [1.0] * len(self.gru)
        self.attrs['slopeTypeIndex'] = xr.DataArray(slope_type_index_values, dims=self.dims['hru'])

        # 3. `hruId` and `gruId` values
        # Note:
        #     If needed, the hard-coded names here can be transformed to be read
        #     from the input objects.
        self.attrs['gruId'] = xr.DataArray(self.gru['COMID'], coords={'gru': self.gru['COMID'].values})
        self.attrs['hruId'] = xr.DataArray(self.hru['COMID'], coords={'hru': self.hru['COMID'].values})

        # 4. `hru2gruId` values
        #    


    def __repr__(
        self,
    ) -> str:
        """Official string representation
        
        FIXME: Needs improvement, for now a sandbox style
        """
        # FIXME: In an ideal world, `self.topology` should be defined
        #        as an object and populated across models. In this
        #        work, we will only focus on main object.
        # forcing information
        forcing_str = f"Forcing files: {self._forcing}"

        # FIXME: In an ideal world, `self.topology` should be defined
        #        as an object and populated. In this work, we will
        #        only focus on main object. Previous assumption was that
        #        `hydrant` can provide such an object.
        # topology information
        # if single site setup is desired, `riv` object
        # will not be necessary; if present, show the number of
        # elements as part of the __repr__
        if self.riv is not None:
            riv_count = len(self.riv)
            riv_str = f"Rivers: {riv_count} "
        else:
            riv_str = f"Rivers: no river network"
        # consider the numbers of `topology['cat']` and
        # `topology['hru']` elements as their __repr__
        cat_count = len(self.cat)
        hru_count = len(self.hru)
        # build strings
        cat_str = f"GRUs: {cat_count}"
        hru_str = f"HRUs: {hru_count}"
        topology_str = cat_str + '\n' + hru_str + '\n' + riv_str

        # object's status
        status = f"Initialized: {self.init}"

        # final string representation is a summary of
        # forcing and topology
        repr_str = forcing_str + '\n' + topology_str + '\n' + status

        return repr_str


