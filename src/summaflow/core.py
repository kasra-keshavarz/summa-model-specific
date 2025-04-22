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
import jinja2 as jij

import pint_xarray
import pint
import pytz

from typing import (
    Dict,
    Sequence,
    Union,
    Any,
    Self,
    Optional,
    Tuple,
    TypeAlias,
)

# built-in libraries
import re
import json
import sys
import glob
import os
import shutil
import warnings

# internal package imports
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
    _calculate_centroids,
    _calculate_polygon_areas,
    _freq_longname,
)
from .geospatial import (
    GeoLayer,
    Stats,
)

# Type definitions
FIDType: TypeAlias = Union[str, int, float]

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
        forcing_name_mapping: Dict[str, str],
        forcing_unit_mapping: Dict[str, str],
        forcing_to_unit_mapping: Dict[str, str],
        topology_data: Dict[str, str | int],
        topology_attrs: Dict[str, str],
        topology_unit_mapping: Dict[str, str],
        topology_to_unit_mapping: Dict[str, str],
        geospatial_data: Optional[Dict[str, Dict]] = None,
        dims: Optional[Dict[str, str]] = default_dims,
        settings: Dict[str, Any] = {},
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
        topology_units : :obj:`dict`
            Original units for all data labels included in `topology_data`
            elements
        topology_to_units : :obj:`dict`
            **Acceptable** units of all data labels included in
            `topology_data` elements SUMMA requires
        topology_attrs : :obj:`dict` of str
            Containing necessary attributes and metadata. **Necessary** keys
            are 
        geospatial_data : :obj:`dict`
            Additional geospatial data used in parameterizing SUMMA.
            Acceptable keys are:
               * *heightCanopyTop*
               * *vegTypeIndex*
               * *soilTypeIndex*
               * *elevation*
        dims : :obj:`dict`, optional
            Critical dimnesion names. Acceptable keys are *hru* and *gru*.

        Note
        ----
        The following lists the elements of `topology_data`
        topology_vars : :obj:`dict` of :obj:`str` keys and values
            Topology variables needed from the objects read through
            `topology_data`. Critical keys are: *segId*, *downSegId*,
            *slope*, *length*, *hruToSegId*, *tan_slope*. Any other
            additional labels will be added to SUMMA's Attribute NetCDF
            file for user's convenience.
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
        self.forcing_vars = forcing_name_mapping
        self._forcing_attrs = forcing_attrs
        self.forcing_units = forcing_unit_mapping
        self.forcing_to_units = forcing_to_unit_mapping
 
        # topology data
        # `riv` is optional, as single-site configurations, do not have
        # any river systems involved
        try:
            self.riv = topology_data.get('riv')
        except:
            self.riv = None

        # `cat` will be added as a virtual property
        try:
            # `gru` can also be called `cat`
            self.gru = topology_data.get('cat')
        except:
            self.gru = topology_data.get('gru')

        # `hru` object
        self.hru = topology_data.get('hru')

        # `topology_*` object
        self.topology_attrs = topology_attrs
        self.topology_units = topology_unit_mapping
        self.topology_to_units = topology_to_unit_mapping

        # Dimension names
        self.dims = dims

        # `init` object to add lazy-like behaviour
        self.init = False

        # Geospatial data
        self.geospatial_data = geospatial_data

        # Workflow settings
        self.settings = settings

        # Pint unit registry
        self._ureg = pint.UnitRegistry(force_ndarray_like=True)

        return

    # custom constructors
    @classmethod
    def from_maf(
        cls,
        layers: Dict[str, Dict] = None,
    ) -> Self:

        return

    # virtual properties
    @property
    def cat(
        self,
    ) -> gpd.GeoDataFrame:
        """aka `gru`"""
        return self.gru

    # special methods
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

    # object methods
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
            dims=default_dims)

        # populating with necessary information obtained upon
        # instantiation; This class variable is a combination of
        # information of `forcing_*` and `topology_*` objects
        #
        # 1. measurement height value through self.forcing_attrs
        # Notes:
        #    1. The `measurement_height` is the assumption of this workflow
        #    that is clearly instructed in the tutorial and API reference.
        #    2. `mHeight` is defined for each `hru` in SUMMA, therefore,
        #       hard-coded here.
        _mheight_name = 'mHeight'
        self.attrs[_mheight_name] = SUMMAWorkflow._mheight(
            forcing_attrs=self._forcing_attrs,
            elements=self.hru,
            element_name=self.dims['hru'],
            height_name='measurement_height',
            height_unit='measurement_height_unit')

        # 2. `slopeTypeIndex` which is "a legacy that is no longer used"
        #     reference: github.com/CH-Earth/CWARHM: step 5/SUMMA/1/1
        # Note:
        #     If needed, the hard-coded names & values here can be transformed
        #     to be read from the input objects.
        _slope_type_index_name = 'slopeTypeIndex'
        self.attrs[_slope_type_index_name] = SUMMAWorkflow._slope_type_index(
            slope_value=int(1),
            elements=self.hru,
            element_name=self.dims['hru'])

        # 3. `hruId` and `gruId` values
        # Note:
        #     If needed, the hard-coded names here can be transformed to be read
        #     from the input objects.
        _gru_id_name = 'gruId'
        _hru_id_name = 'hruId'
        self.attrs[_gru_id_name] = SUMMAWorkflow._elements(
            self.gru,
            self.topology_attrs['gru_fid'],
            'gru',
            self.gru[gru_fid])
        self.attrs[_hru_id_name] = SUMMAWorkflow._elements(
            self.hru,
            hru_fid,
            'hru',
            self.hru[hru_fid])

        # 4. `hru2gruId` values
        # Note:
        #     The `topology` needs to become a systematic object using
        #     `hydrant` but unfortunately the project does not have much
        #     support.
        _hru_mapping_gru_name = 'hru2gruId'
        mapping = SUMMAWorkflow._mapping_hru(
            gru=self.gru,
            hru=self.hru,
            gru_fid=gru_fid,
            hru_fid=hru_fid,
            mapping_fid=gru_fid,
        )
        self.attrs[_hru_mapping_gru_name] = xr.DataArray(
            pd.Series(mapping),
            dims=["hru"],    # Name of the dimension
            name="hru2gruId" # Name of the variable
        )

        # 5. `latitude` and `longitude`
                # Note:
        #     The values are extracted for the centroid of each `hru`
        # First, calculate centroids; EPSG 6933 is hard-coded for accuracy
        # the returned object will contain data labels `centroid_y` and
        # `centroid_x`
        coords_defs = {
            'latitude': 'centroid_x',
            'longitude': 'centroid_y',
        }
        centroids = _calculate_centroids(self.hru)
        for k, v in coords_defs.items():
            self.attrs[k] = xr.DataArray(
                centroids[v],
                coords={'hru': centroids[hru_fid]})

        # 6. `area` values
        # Notes:
        #    the `target_area_unit` is hard-coded as it is SUMMA's
        #    default requirement. The returned object of this function
        #    includes a column named `area`.
        #    The unit of area is critical, so they have been assigned now;
        #    for all other variables, they are added later.
        areas = _calculate_polygon_areas(self.hru, target_area_unit='m^2')
        area_unit = areas['area'].pint.units
        self.attrs['HRUarea'] = xr.DataArray(
            data=areas['area'].pint.magnitude,
            coords={'hru': areas[hru_fid].values},
            attrs={'unit': area_unit})

        # Up to Martyn & his staff: the values decided in various workflows
        #    are not similar.
        #
        # 7. `tan_slope` values
        #    Darri: uses np.gradient
        #    Mohamed: uses river slope values
        #    Wouter: uses 0.1 constant value
        #
        # 8. `contourLength` values
        #    Darri: uses the length of intersection between an HRU
        #           and its downstream neighbor. [KK: how one can calculate
        #           a downstream HRU?]. Or, he uses the square root of HRU
        #           area.
        #    Mohamed: the length around the riparian zone (i.e., length of
        #             stream *2).
        #    Wouter: uses 30m constant value
        #
        # 9. `downHRUindex` values
        #    Martyn: for non-contiguous HRUs, assign downHRUindex the
        #            riparian, contiguous one
        #    Darri: sets the downHRUindex based on elevation data
        #    Mohamed: all constant 0 values
        #    Wouter: Set the downHRUindex based on elevation data

        # 10. `geospatial` layers
        # 10.1 `eleveation` layer
        _elv_name = 'elevation'
        elv = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_elv_name],
            _elv_name,
            'mean',
            'hru')
        # 10.2 `vegTypeIndex` layer
        _veg_name = 'vegTypeIndex'
        veg = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_veg_name],
            _veg_name,
            'majority',
            'hru')
        # 10.3 `soilTypeIndex` layer
        _soil_name = 'soilTypeIndex'
        soil = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_soil_name],
            _soil_name,
            'majority',
            'hru')

        _geolayers = {
            _elv_name: elv,
            _veg_name: veg,
            _soil_name: soil
        }
        
        # Adding Geospatial layers to `self.attrs`
        self.attrs.update(_geolayers)

        # Updating local attributes 
        # Assign attributes to each variable
        for var_name, attrs in attributes_local_attrs_default.items():
            if var_name in self.attrs:
                self.attrs[var_name].attrs.update(attrs)
        # Updating global attributes
        self.attrs = self.attrs.assign_attrs(
            attributes_global_attrs_default
        )

        return self.attrs

    def init_forcing(
        self,
    ) -> None:
        """Prepare forcing dataset for the SUMMA setup. The preparation step
        involves name change, unit adjustments, and sorting (SUMMA-specific)
        NetCDF files as the forcing inputs. The default inputs are always
        NetCDF files, so no flexibility is implemented for other data formats.
        Another step is to sort the NetCDF files and create a textual file
        listing all inputs.

        Notes
        -----
        - Timezone naming scheme (TZ) follows IANA's convention found at the
          following link (version 2025b):
            https://data.iana.org/time-zones/releases/tzdb-2025b.tar.lz
 
        """
        # Iterate over the dataset files and implement necessary adjustments
        # FIXME: Can be parallelized via Dask, though creating dask clusters
        #        on HPCs can be challenging.
        if self._forcing_attrs['forcing_time_zone']:
            _forcing_tz = self._forcing_attrs['forcing_time_zone'].lower()
        else:
            _forzing_tz = 'local'

        if self._forcing_attrs['target_time_zone']:
            _target_tz = self._forcing_attrs['target_time_zone'].lower()
        else:
            _target_tz = 'local'

        # Specify the `forcing` and `target` timezones and also the string
        # for the "fileManager"
        forcing_tz, target_tz, tz_info = SUMMAWorkflow._specify_tz(_forcing_tz, _target_tz)
        self._forcing_attrs['tz_info'] = tz_info

        # Time travel to late-90s and go over the NetCDF files one
        # by one

        # FIXME: Considering Dask for parallelization in near future
        for forcing in self._forcing:
            # First read the forcing file, using Xarray
            ds = xr.open_dataset(forcing)

            # Change timezone and assing tz-naive datetime64[ns] timestamps
            if target_tz not in ('local'):
                ds = ds.assign_coords({
                   'time': ds.time.to_index().tz_localize(forcing_tz).tz_convert(target_tz).tz_localize(None)
                })

            # Rename variables
            ds = ds.rename_vars(self.forcing_vars)

            # Change units
            # Check if all units provided are present in the forcing file(s)
            for k in self.forcing_units:
                if k not in ds:
                    raise ValueError(f"item {k} defined in "
                                     "`forcing_unit_mapping` cannot be found"
                                     " in `forcing_name_mapping` values.")

            # If all elements of `variables` not found in `units`,
            # assign them to None
            for v in ds:
                if v not in self.forcing_vars.values():
                    self.forcing_units[v] = 'dimensionless'

            # Now assign the units
            ds = ds.pint.quantify(units=self.forcing_units, unit_registry=self._ureg)
 
            # If `to_units` is defined
            if self.forcing_to_units:
                ds = ds.pint.to(units=self.forcing_to_units)

            # Print the netCDF file
            ds = ds.pint.dequantify()

        return

    def _unit_change(
        self: Self,
    ) -> xr.Dataset:

        return

    # static methods
    @staticmethod
    def _specify_tz(
        forcing: str,
        target: str,
    ) -> Tuple[str, str, str]:
        """Based on SUMMA's functionality, return a tuple containing
        `forcing` and `target` time-zones that can be used with pandas.Index
        `.tz_localize(...)` and `.tz_convert(...)` functionality.

        Parameters
        ----------
        forcing : str
            User-specified reference, original forcing dataset timezone. Can
            be ``UTC``, ``GMT``, ``local`` or an IANA timezone string.
        target : str
            Target timezone for SUMMA configurations. Can be ``UTC``, ``GMT``,
            ``local`` or an IANA timezone string.

        Returns
        -------
        Tuple
            A tuple object containing forcing and target timezones for further
            adjustments by `.tz_localize(...)` and `.tz_convert(...)`
            functions. The last element is `tz_info` that should be specified
            in SUMMA's "filemanager".

        Raises
        ------
        ValueError
            - If forcing is set to `local` but the target is another known
              timezone, it is considered ambiguous.
            - If target timezone is ambiguous.

        Notes
        -----
        - If `forcing` and `target` cannot be interpretted, or not acceptable
          by the model, necessary raises are thrown.
        """
        # Routine error checks
        if not isinstance(forcing, str):
            raise TypeError("forcing time zone needs to be of dtype `str`.")
        if not isinstance(target, str):
            raise TypeError("target time zone needs to be of dtype `str`.")

        # Working with lowercase strings
        target = target.lower()
        forcing = forcing.lower()

        # List of lowercase IANA timezones
        time_zones = (tz.lower() for tz in pytz.common_timezones)

        # `local` timezone is ambiguous alone, unless both `forcing` and
        # `target` are set to `local`, trusting user's time manipulations.
        if forcing in ("local") and target not in ("local"):
            raise ValueError("If forcing time zone is set to `local`, the" 
                " target must necessarily be `local`.")

        # If set to `utc` (or `gmt`), go and adjust `forcing`
        # if necessary
        if target in ('utc', 'gmt'):
            # If target and forcing tzs are both UTC, do nothing
            if forcing in ('utc', 'gmt'):
                tz_info = 'utcTime'

            # If target is `utc` and forcing is a local tz, convert to
            # `utc`
            elif forcing in time_zones:
                tz_info = 'utcTime'

            else:
                raise ValueError("forcing time zone ambiguous.")

        # If user knows what a target local timezone is and can handle it
        # no need to change the tzs
        elif target in ('local'):
            if forcing not in ('local'):
                raise ValueError("target time zone ambiguous.")
            else:
                tz_info = 'localTime'

        # If target timezone is set to a IANA-standard timezone
        elif target in time_zones:
            if forcing in time_zones or forcing in ('utc' , 'gmt'):
                # If user specifies an IANA-standard timezone
                tz_info = 'localTime'
            else:
                raise ValueError("forcing time zone ambiguous.")

        # Else, if user provided invalid value for target time-zone,
        # raise ValueError
        else:
            raise ValueError("target time zone ambiguous.")

        return (forcing, target, tz_info)

    @staticmethod 
    def _specify_time_encodings(
        ds: xr.Dataset,
        time_var: str = 'time',
    ) -> Dict[str, Dict[str, str]]:
        """Necessary adhoc modifications on the forcing object's
        encoding
        """
        # empty encoding dictionary of the `time` variable
        ds[time_var].encoding = {}
        # estimate the frequency offset value
        _freq = pd.infer_freq(ds[time_var])
        # get the full name
        _freq_long = utils.freq_long_name(_freq)

        # Encoding dictionary appropriate to be included as a local attribute
        _encoding = {
            'time': {
                'units': f'{_freq_long} since 1900-01-01 12:00:00'
            }
        }

        return _encoding

    @staticmethod
    def _mheight(
        forcing_attrs: Dict[str, str],
        elements: Sequence[str | float | int],
        element_name: str = 'hruId',
        height_name: str = 'measurement_height',
        height_unit: str = 'measurement_height_unit',
    ) -> xr.DataArray:
        """Returning SUMMA-specific *mHeight* variable based on [1].

        Parameters
        ----------
        forcing_attrs : Dict[str, str]
            A Dictionary describing necessary attribute of forcing variables.
        elements : Sequence of str, float, or int
            A seuqence of GRU values.
        element_name : str
            Name of the `elements` being used as the dimension of returned 
            |DataArray|.
        height_name : str
            Variable defining measurement height value in `forcing_attrs` keys.

        References
        ----------
        .. [1] https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/#local-attributes-file
        """
        if not isinstance(forcing_attrs, dict):
            raise TypeError("`forcing_attrs` must be a dictionary of string values.")

        mHeight_name = forcing_attrs.get(height_name)
        mHeight_values = [forcing_attrs.get(height_name)] * len(elements)
        attrs = {'unit': forcing_attrs.get(height_unit)}

        return xr.DataArray(mHeight_values, dims=element_name, attrs=attrs)

    @staticmethod
    def _slope_type_index(
        slope_value : int,
        elements : Sequence[int | str | float],
        element_name : str = 'hruId',
    ) -> xr.DataArray:
        """Returning SUMMA-specific *slopeTypeIndex* variable based on [1].
        Based on [1], this value is a legacy variable and can be set to an
        integer value of 1 and repeated for all computational elements.

        Parameters
        ----------
        slope_value : Dict[str, str]
            A Dictionary describing necessary attribute of forcing variables.
        elements : Sequence of str, float, or int
            A seuqence of HRU/GRU values.
        element_name : str
            Name of the `elements` being used as the dimension of returned 
            |DataArray|.

        References
        ----------
        .. [1] https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/#local-attributes-file
        """
        if not isinstance(slope_value, int | float):
            raise TypeError("`slope_value` must be an integer or a float")
        slope_type_index_values = [int(slope_value)] * len(elements)

        return xr.DataArray(slope_type_index_values, dims=element_name)

    @staticmethod
    def _tan_slope(
        hru: gpd.GeoDataFrame,
    ) -> xr.DataArray:
        """To be completed after discussion"""

        return

    @staticmethod
    def _contour_length(
        hru: gpd.GeoDataFrame,
    ) -> xr.DataArray:
        """To be completed after discussion"""

        return

    @staticmethod
    def _down_hru_index(
        hru: gpd.GeoDataFrame,
    ) -> xr.DataArray:
        """To be completed after discussion"""

        return

    @staticmethod
    def _elements(
        geom: gpd.GeoDataFrame,
        fid: str,
        dim_name: str,
        dim_value: str,
    ) -> xr.DataArray:
        """Return element values in from of a |DataArray|

        Parameters
        ----------

        Returns
        -------
        """
        if not isinstance(geom, gpd.GeoDataFrame):
            raise ValueError("`geom` must be a geopandas.GeoDataFrame.")

        return xr.DataArray(geom[fid], coords={dim_name: dim_value})

    @staticmethod
    def _mapping_hru(
        gru: gpd.GeoDataFrame,
        hru: gpd.GeoDataFrame,
        gru_fid: str,
        hru_fid: str,
        mapping_fid: str,
    ) -> Dict[FIDType, FIDType]:
        """Create a mapping dictionary from HRU identifiers to their corresponding GRU identifiers.

        Parameters
        ----------
        gru : gpd.GeoDataFrame
            GeoDataFrame containing GRU information. Must contain the columns 
            specified by `gru_fid`.
        hru : gpd.GeoDataFrame
            GeoDataFrame containing HRU information. Must contain the columns
            specified by `hru_fid`.
        gru_fid : str
            Column name in `gru` that contains the GRU (Grouped Response Unit) identifiers.
        hru_fid : str
            Column name in `hru` that contains the HRU (Hydrologic Response Unit) identifiers.
        mapping_fid : dict of FIDType keys and values
            Name mapping information from `hru` to `gru` data within `hru` object.

        Returns
        -------
        Dist[FIDType, FIDType]
            A dictionary where keys are HRU identifiers and values are the corresponding
            GRU identifier for each HRU.

        Raises
        ------
        TypeError
            If input hru is not a GeoDataFrame
            If gru_fid or hru_fid are not strings
        ValueError
            If specified columns are not in the GeoDataFrame
            If there are duplicate HRU identifiers
            If any GRU values are null/missing

        Notes
        -----
        - The HRU and GRU defintions can be found in [1].

        References
        ----------
        .. [1] Clark et al., (2015). “A Unified Approach for Process-Based
               Hydrologic Modeling: 1. Modelling Concept”
               DOI: 10.1002/2015WR017198
        """
        # Input type validation
        if not isinstance(hru, gpd.GeoDataFrame):
            raise TypeError("`hru` must be a geopandas.GeoDataFrame object.")
        if not isinstance(gru, gpd.GeoDataFrame):
            raise TypeError("`gru` must be a geopandas.GeoDataFrame object.")
        if not isinstance(gru_fid, str) or not isinstance(hru_fid, str):
            raise TypeError("`gru_fid` and `hru_fid` must be of type string.")

        # Column existence validation
        missing_cols = [
            f"{hru_fid} (in hru)" if hru_fid not in hru.columns else None,
            f"{gru_fid} (in gru)" if gru_fid not in gru.columns else None
        ]
        missing_cols = [col for col in missing_cols if col is not None]
        if missing_cols:
            raise ValueError(f"Missing columns: {', '.join(missing_cols)}")

        # Check to see if `hru` and `gru` are equal objects
        if hru.equals(gru):
            pass
        else:
            # Data quality checks for HRU and GRU identifiers
            for gdf, fid, name in [(hru, hru_fid, "HRU"), (gru, gru_fid, "GRU")]:
                if gdf[fid].duplicated().any():
                    raise ValueError(f"Duplicate {name} identifiers found - each {name} should be unique.")
                if gdf[fid].isnull().any():
                    raise ValueError(f"Null values found in {name} identifiers.")

        # The values associated with mapping_fid should have equal length compared
        # to the corresponding `gru_fid` values
        _gru_elements = set(gru[gru_fid])
        _mapping_elements = set(hru[mapping_fid])
        if not _gru_elements.issubset(_mapping_elements):
            raise IndexError("`hru` does not cover all available `gru` elements.")

        # Create and return the mapping
        mapping = hru.groupby(hru_fid)[mapping_fid].first().to_dict()
 
        # Verify mapping completeness
        if len(mapping) != len(hru):
            raise RuntimeError("Unexpected error in mapping creation - size mismatch")
 
        return mapping

    @staticmethod
    def _geolayer_info(
        layer: GeoLayer,
        layer_name: str,
        stat_names: Sequence[str] | str,
        dim_name: str,
        unit: pint.Unit = None,
    ) -> xr.DataArray:
        """GeoLayer information in form of a |DataArray|

        Parameters
        ----------
        layer : GeoLayer
            A GeoLayer holding necessary information needed for a SUMMA setup.
        name : str
            The name of the GeoLayer assigned to the returned |DataArray|
            needed for a SUMMA setup.
        stat_names : Sequence of str or str
            One or more statistics that needs to be reported for a SUMMA
            setup.
        dim_name : str
            The dimension name used for representing values for each SUMMA
            element.
        unit : pint.Unit or `None`
            If GeoLayer's unit needs to be changed, this sepcifies the target
            unit of interest.

        Returns
        -------
        |DataArray|
            A |DataArray| describing the layer of interest needed for a SUMMA
            setup.

        Notes
        -----
        - It is assumed that layer holds relevant unit and statistics,
          however, upon need, that can change into `unit`.
        """
        # Routine checks
        if not isinstance(layer, GeoLayer):
            raise TypeError("`layer` must be of type GeoLayer.")

        # If `unit` is defined, convert the unit of the layer
        if unit:
            layer.to_unit(unit)

        da = xr.DataArray(
            data=layer.stats[stat_names],
            coords={
                dim_name: layer.stats[stat_names].index.values,
            },
            attrs={'unit': layer.unit},
            name=layer_name,
        )

        return da
