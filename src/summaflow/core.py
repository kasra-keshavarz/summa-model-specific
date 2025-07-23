"""
This package automates SUMMA model setup through a flexible workflow that can
be set up using a JSON configuration file, a Command Line Interface (CLI),
directly inside a Python script/environment, or other commonly used
configuration interfaces
"""

# Built-in libraries
import itertools
import math
import os
import pathlib
import shutil
import warnings
import json
import glob
import re
import numbers

from importlib.resources import files, as_file

from typing import (
    Any,
    Dict,
    Optional,
    Sequence,
    Tuple,
    TypeAlias,
    Union,
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
import pytz
import xarray as xr
import pint_xarray

from jinja2 import Environment, FileSystemLoader

# Internal package imports
from . import utils  # ./utils.py file
from ._default_dicts import (
    attributes_global_attrs_default,
    attributes_local_attrs_default,
    cold_state_global_attrs_default,
    cold_state_local_attrs_default,
    default_dims,
    forcing_global_attrs_default,
    forcing_local_attrs_default,
    trial_params_local_attrs_default,
    trial_params_global_attrs_default,
    model_decisions_default,
)
from .geospatial import (
    GeoLayer,
    Stats,
)
from .__about__ import __version__
from .logging_config import setup_logger

# Type definitions
# Path type
try:
    from os import PathLike
except ImportError:  # for Python < 3.6
    PathLike = str

# GIS object FID type
FIDType: TypeAlias = Union[str, int, float]

# Module-level logger
logger = setup_logger(__name__, level=10)


class SUMMAWorkflow:
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
        forcing_data: Sequence[str | os.PathLike] = [],
        forcing_attrs: Dict[str, str] = {},
        forcing_name_mapping: Dict[str, str] = {},
        forcing_unit_mapping: Dict[str, str] = {},
        forcing_to_unit_mapping: Dict[str, str] = {},
        topology_data: Dict[str, str | int] = {},
        topology_attrs: Dict[str, str] = {},
        topology_unit_mapping: Dict[str, str] = {},
        topology_to_unit_mapping: Dict[str, str] = {},
        cold_state: Dict[str, Any] = {},
        geospatial_data: Optional[Dict[str, Dict]] = None,
        settings: Dict[str, Any] = {},
        decisions: Dict[str, str] = {},
        fillna: Dict[str, Dict] = {},
        dims: Optional[Dict[str, str]] = default_dims,
        **kwargs,
    ) -> None:
        """
        Initialize a SUMMA workflow configuration. All input arguments are
        considered optional to ease independence of various configuration
        stages.

        Parameters
        ----------
        forcing_data : Sequence[str | os.PathLike], optional
            Sequence of paths to NetCDF files containing forcing data for SUMMA.
        forcing_attrs : Dict[str, str], optional
            Dictionary of attributes for forcing variables. Must include
            'measurement_height' and 'measurement_height_unit' keys.
        forcing_name_mapping : Dict[str, str], optional
            Mapping between variable names in forcing files and SUMMA-accepted names.
            Keys are source names, values are SUMMA-accepted names.
        forcing_unit_mapping : Dict[str, str], optional
            Original units for forcing variables. Keys are SUMMA-accepted names,
            values are original units (following Pint's default_en.txt standards).
        forcing_to_unit_mapping : Dict[str, str], optional
            Target units for forcing variables. Keys are SUMMA-accepted names,
            values are required SUMMA units (following Pint's default_en.txt standards).
        topology_data : Dict[str, str | int], optional
            Topology data including river, catchment, and HRU information.
            Expected keys:
            - 'riv': Path or GeoDataFrame for river network
            - 'cat': Path or GeoDataFrame for catchment polygons (GRUs)
            - 'hru': Path or GeoDataFrame for HRU polygons
        topology_attrs : Dict[str, str], optional
            Attributes for topology data. Must include:
            - 'gru_fid': Common label between river and catchment data
            - 'hru_fid': Mapping between HRU and GRU values
        topology_unit_mapping : Dict[str, str], optional
            Original units for topology variables.
        topology_to_unit_mapping : Dict[str, str], optional
            Target units for topology variables required by SUMMA.
        cold_state : Dict[str, Any], optional
            Dictionary defining initial conditions for model state variables.
            Expected structure:

            .. code-block:: python

                {
                    'layers': {
                        'nSoil': int,  # Number of soil layers (e.g., 8)
                        'nSnow': int   # Number of snow layers (e.g., 0 for no snow)
                    },
                    'states': {
                        # Scalar states (single global value or sequence of values per HRU)
                        'scalarCanopyIce': float       # Initial canopy ice [kg/m²]
                        'scalarCanopyLiq': float       # Initial canopy liquid [kg/m²]
                        'scalarSnowDepth': float       # Initial snow depth [m]
                        'scalarSWE': float             # Initial snow water equivalent [m]
                        'scalarSfcMeltPond': float     # Surface melt pond storage [m]
                        'scalarAquiferStorage': float  # Initial aquifer storage [m]
                        'scalarSnowAlbedo': float      # Initial snow albedo [0-1]
                        'scalarCanairTemp': float      # Canopy air temperature [K]
                        'scalarCanopyTemp': float      # Canopy temperature [K]

                        # Layer states (array per layer and HRU)
                        'mLayerTemp': float | List[float]       # Layer temperatures [K]
                        'mLayerVolFracIce': float | List[float] # Volumetric ice fraction [m³/m³]
                        'mLayerVolFracLiq': float | List[float] # Volumetric liquid fraction [m³/m³]
                        'mLayerMatricHead': float | List[float] # Matric head [m]
                        'mLayerDepth': List[float]              # Layer depths [m] (e.g., [0.025, 0.075,...])
                    }
                }

            Typical cold start values assume:
            - No snow cover (nSnow=0, scalarSWE=0)
            - Moist soil (mLayerVolFracLiq=0.4)
            - 10°C initial temperature (283.16K)
            - Moderate aquifer storage (0.4m)
        geospatial_data : Optional[Dict[str, Dict]], optional
            Geospatial attributes for SUMMA parameterization. Expected structure:

            .. code-block:: python

                {
                    'elevation': GeoLayer(**{
                        'stats': Stats(ArrayLike),  # Must include `mean` index
                        'unit': str
                    }),
                    'soilTypeIndex': {
                        'stats': Stats(ArrayLike),  # Must include `majority` index
                        'unit': str
                    }),
                    'vegTypeIndex': {
                        'stats': Stats(ArrayLike),  # Must include `majority` index
                        'unit': str
                    }),
                    'tan_slope': {
                        'stats': Stats(ArrayLike),  # Must include `mean` index
                        'unit': str
                    }),
                    'contourLength': {
                        'stats': Stats(ArrayLike),  # Must include `length` index
                        'unit': str
                    }),
                    'downHRUindex': {
                        'stats': Stats(ArrayLike),  # Must include `index` index
                        'unit': str
                    }),
                }

        settings : Dict[str, Any], optional
            Model configuration settings. Structure:

            .. code-block:: python

                {
                    'model_path': str,           # Path to SUMMA settings directory
                    'start_date': str,           # Simulation start datetime
                    'end_date': str,             # Simulation end datetime
                    'verbose': bool,             # Enable verbose output
                }

        decisions : Dict[str, str], optional
            Model process decisions. The keys are SUMMA-accepted model options
            and the corresponding values are available modelling decisions. If
            a decision is not provided, default SUMMA values will be
            consdiered.
        fillna : Dict[str, Dict], optional
            Fill values for missing data.

            .. code-block:: python

                {
                    'geospatial_data': {
                        'elevation': float,    # Fill value for missing elevation
                        'soilTypeIndex': int,  # Fill value for missing soil types
                        'vegTypeIndex': int,   # Fill value for missing vegetation
                        'tan_slope': float,    # Fill value for missing slopes
                        'contourLength': float # Fill value for missing contours
                    }
                }

        dims : Optional[Dict[str, str]], optional
            Critical dimension names. Expected keys:
            - 'hru': Hydrologic Response Unit dimension
            - 'gru': Grouped Response Unit dimension
        **kwargs
            Additional keyword arguments.

        Notes
        -----
        Topology Data Requirements:
        - 'topology_vars' dictionary should include these critical variables:
            * 'segId': Segment identifiers
            * 'downSegId': Downstream segment identifiers
            * 'slope': Segment slope
            * 'length': Segment length
            * 'hruToSegId': HRU to segment mapping
            * 'tan_slope': Tangent of slope
        - If 'hru' is not provided, each HRU will be equivalent to GRU elements.
        Geospatial Data Requirements:
        - ``Stats`` object requirements:
            * elevation: must include ``mean`` stat
            * soilTypeIndex: must include ``majority`` stat
            * vegTypeIndex: must include ``majority`` stat
            * tan_slope: must include ``mean`` stat
            * contourLength: must include ``length`` stat
            * downHRUindex: must include ``index`` stat
        Unit Handling:
        - All unit specifications should follow Pint's default_en.txt standards.
        - Unit conversion is automatically handled when both original and target units are provided.
        - In cases where target unit cannot be specified, this conversion will
          be taken care of internally.
        """
        # assign necessary attributes
        # FIXME: This needs to turn into its own object, but for the
        #        sake of timing of this deliverable, we compromise and treat
        #        them as different variables with `forcing_` prefix.
        # `forcing_data` needs to be a sequence of paths
        if isinstance(forcing_data, Sequence) and not \
            isinstance(forcing_data, (bytes, str)):
            self._forcing = sorted(forcing_data)
        else:
            self._forcing = None

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

        # Cold State (initial condition) data
        self.cold_state_data = cold_state

        # Geospatial data
        self.geospatial_data = geospatial_data

        # Workflow settings
        self.settings = settings
        # Make the directory for the outputs
        if 'model_path' in self.settings:
            try:
                os.makedirs(self.settings['model_path'], exist_ok=True)
            except ValueError as e:
                raise ValueError(f"Invalid model path: {self.settings['model_path']}") from e
            except OSError as e:
                pass

        # Pint unit registry
        self._ureg = pint.UnitRegistry(force_ndarray_like=True)

        # Auxillary data dictionary for intermediate important information
        if 'auxillary' in kwargs.keys():
            self.auxillary = kwargs['auxillary']
        else:
            self.auxillary = {}

        # `init` object to add lazy-like behaviour
        self.auxillary['init'] = []

        # If `dt_init` not in `auxillary`, warn defaulting to forcing time-step
        if 'dt_init' not in self.auxillary:
            warnings.warn("`dt_init` not provided in auxillary dictionary;"
                " defaulting to forcing time-step.")

        # Jinja2 environment for templating
        def raise_helper(msg):
            raise Exception(msg)

        self.package_root = files('summaflow')

        self.jinja2_env = Environment(
            loader=FileSystemLoader(self.package_root.joinpath('templates')),
            trim_blocks=True,
            lstrip_blocks=True,
            line_comment_prefix='##',
        )
        self.jinja2_env.globals['raise'] = raise_helper

        # Registering model decisions
        if decisions:
            self.decisions = decisions
        else:
            self.decisions = None

        # Registering fillna
        self.fillna = fillna

        # Defining verbosity variable
        self.verbose = self.settings.get('verbose', False)
        if self.verbose:
            logger.info("SUMMA workflow initialized")

    # custom constructors
    @classmethod
    def from_maf(
        cls,
        forcing_data: PathLike | str = None, # type: ignore
        topology_data: Dict[str, PathLike | str] = {}, # type: ignore
        geospatial_data: Optional[Dict[str, Dict]] = None, # type: ignore
        **kwargs,
    ) -> Self:
        """Construct a SUMMAWorkflow object from MAF-compatible files"""
        logger.info("Constructing SUMMAWorkflow from MAF files")

        forcing_to_unit_mapping = {
            'pptrate': 'millimeter / second',
            'airtemp': 'kelvin',
            'airpres': 'pascal',
            'LWRadAtm': 'watt / meter ** 2',
            'SWRadAtm': 'watt / meter ** 2',
            'spechum': 'dimensionless',
            'windspd': 'meter / second',
        }

        # building forcing data list
        forcing_data_obj = glob.glob(os.path.join(forcing_data, '**/*.nc*'), recursive=True)

        # topology data
        topology_data_obj = {k: gpd.read_file(v) for k, v in topology_data.items() if k in ['riv', 'cat', 'hru']}

        # geospatial data
        # layers needed by the setup workflow
        # elevation
        elv = GeoLayer.from_maf(
            maf_stats=geospatial_data['elevation'],
            maf_layer=None,
            maf_geolayer=os.path.join(topology_data['cat']),
            unit = 'meters',
        )
        # landcover
        landcover = GeoLayer.from_maf(
            maf_stats=geospatial_data['vegTypeIndex'],
            maf_layer=None,
            maf_geolayer=os.path.join(topology_data['cat']),
            unit = 'dimensionless',
        )
        # USDA soil classes
        soil = GeoLayer.from_maf(
            maf_stats=geospatial_data['soilTypeIndex'],
            maf_layer=None,
            maf_geolayer=os.path.join(topology_data['cat']),
            unit = 'dimensionless',
        )

        # custom layers for `tan_slope`, `contourLength` and `downHRUindex`
        # until relevant workflows are implemented inside `gistool`--sorry
        # For now, look at various constructors for "GeoLayer"
        slope = GeoLayer( # workflow needs `mean` stat
            stats=Stats(pd.DataFrame([0.1] * len(topology_data_obj['cat']), index=topology_data_obj['cat']['COMID'], columns=['mean'])),
            unit='dimensionless',
        )

        contour = GeoLayer( # workflow needs `length` stat
            stats=Stats(
                pd.DataFrame(
                    topology_data_obj['cat'].to_crs('ESRI:54009').length.to_list(), index=topology_data_obj['cat']['COMID'], columns=['length'])),
            unit='meter',
        )

        hru_index = GeoLayer( # workflow needs `index` "stat"
            stats=Stats(pd.DataFrame([0] * len(topology_data_obj['cat']), index=topology_data_obj['cat']['COMID'], columns=['index'])),
            unit='dimensionless',
        )

        # build the geospatial_data_obj dictionary
        geospatial_data_obj = {
            'elevation': elv,
            'vegTypeIndex': landcover,
            'soilTypeIndex': soil,
            'tan_slope': slope,
            'contourLength': contour,
            'downHRUindex': hru_index
        }

        # build the class instance
        instance = cls(
            forcing_data=forcing_data_obj,
            forcing_to_unit_mapping=forcing_to_unit_mapping,
            topology_data=topology_data_obj,
            geospatial_data=geospatial_data_obj,
            **kwargs
        )

        return instance

    @classmethod
    def from_dict(
        cls: 'SUMMAWorkflow',
        init_dict: Dict = {},
    ) -> 'SUMMAWorkflow':
        """
        Constructor to use a dictionary to instantiate
        """
        if len(init_dict) == 0:
            raise KeyError("`init_dict` cannot be empty")
        assert isinstance(init_dict, dict), "`init_dict` must be a `dict`"

        return cls(**init_dict)

    @classmethod
    def from_json(
        cls: 'SUMMAWorkflow',
        json_str: str,
    ) -> 'SUMMAWorkflow':
        """
        Constructor to use a loaded JSON string
        """
        # building customized SUMMAWorkflow's JSON string decoder object
        decoder = json.JSONDecoder(object_hook=SUMMAWorkflow._json_decoder)
        json_dict = decoder.decode(json_str)
        # return class instance
        return cls.from_dict(json_dict)

    @classmethod
    def from_json_file(
        cls: 'SUMMAWorkflow',
        json_file: 'str',
    ) -> 'SUMMAWorkflow':
        """
        Constructor to use a JSON file path
        """
        with open(json_file) as f:
            json_dict = json.load(f, object_hook=SUMMAWorkflow._json_decoder)

    @classmethod
    def _from_maf_json_file(
        cls: 'SUMMAWorkflow',
        json_file: 'str',
    ) -> 'SUMMAWorkflow':
        """
        Constructor to use a MAF-compatible JSON file path
        """
        with open(json_file) as f:
            json_dict = json.load(f, object_hook=SUMMAWorkflow._json_decoder)

        # build the class instance
        return cls.from_maf(**json_dict)

    @staticmethod
    def _env_var_decoder(s):
        """
        OS environmental variable decoder
        """
        # RE patterns
        env_pat = r'\$(.*?)/'
        bef_pat = r'(.*?)\$.*?/?'
        aft_pat = r'\$.*?(/.*)'
        # strings after re matches
        e = re.search(env_pat, s).group(1)
        b = re.search(bef_pat, s).group(1)
        a = re.search(aft_pat, s).group(1)
        # extract environmental variable
        v = os.getenv(e)
        # return full: before+env_var+after
        if v:
            return b+v+a
        return s

    @staticmethod
    def _json_decoder(obj):
        """
        Decoding typical JSON strings returned into valid Python objects
        """
        if obj in ["true", "True", "TRUE"]:
            return True
        elif obj in ["false", "False", "FALSE"]:
            return False
        elif isinstance(obj, str):
            if '$' in obj:
                return SUMMAWorkflow._env_var_decoder(obj)
            if SUMMAWorkflow._is_valid_integer(obj):
                return int(obj)
        elif isinstance(obj, dict):
            return {SUMMAWorkflow._json_decoder(k): SUMMAWorkflow._json_decoder(v) for k, v in obj.items()}
        return obj
 
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
        forcing_files_len = len(self._forcing)
        forcing_str = f"Forcing files: {forcing_files_len} files"

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
            riv_str = f"Rivers segments: {riv_count} "
        else:
            riv_str = "Rivers segments: no river network"
        # consider the numbers of `topology['cat']` and
        # `topology['hru']` elements as their __repr__
        cat_count = len(self.cat)
        hru_count = len(self.hru)
        # build strings
        cat_str = f"GRU count: {cat_count}"
        hru_str = f"HRU count: {hru_count}"
        topology_str = cat_str + '\n' + hru_str + '\n' + riv_str

        # object's status
        if self.auxillary['init']:
            status = f"Initialized: {self.auxillary['init']}"
        else:
            status = "Initialized: False"

        # final string representation is a summary of
        # forcing and topology
        repr_str = forcing_str + '\n' + topology_str + '\n' + status

        return repr_str

    # instance methods
    def init_attrs(
        self,
        return_ds: bool = False,
        save: bool = False,
        save_path: Optional[PathLike | str] = None,
    ) -> None:
        """Initialize the necessary objects for the experient

        Parameters
        ----------
        returns_ds: bool, defaults to `False`
            Returing the |Dataset| if required.
        save: bool, defaults to False
            Saving the resulting |Dataset|.
        save_path: `PathLike` or `str`
            Path to save the resulting |Dataset| including the filename.

        Returns
        -------
        None

        """
        # variables to build an empty SUMMA-specific `attribute` object
        # since `self.topology_attrs['fid']` assumes to be found in both.
        # We use `self.gru` to build the `gru` variable values
        # the keys to `variables` is mandated by `utils._init_empty_ds`
        # function
        if self.verbose:
            logger.info("Initializing attributes for SUMMA workflow...")

        if 'attrs' not in self.auxillary['init']:
            self.auxillary['init'].append('attrs')
        else:
            del self.attrs

        # If necessary fields are not provided, raise an error
        if 'measurement_height' not in self._forcing_attrs:
            raise ValueError("`measurement_height` is a required field in `forcing_attrs`.")
        if 'measurement_height_unit' not in self._forcing_attrs:
            raise ValueError("`measurement_height_unit` is a required field in `forcing_attrs`.")

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

        if self.verbose:
            logger.info("Creating attributes xarray.Dataset")
        # create an empty xarray.Dataset to be populated with SUMMA-specific
        # attributes
        self.attrs = utils._init_empty_ds(
            variables=variables,
            dims=default_dims) # global variable imported from ._default_dicts

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
        if self.verbose:
            logger.info("Adding `mHeight` attribute")
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
        #   - If needed, the hard-coded names & values here can be transformed
        #     to be read from the input objects.
        #   - If provided through the `geospatial_data`, it will be
        #     overwritten.
        if self.verbose:
            logger.info("Adding `slopeTypeIndex` attribute")
        _slope_type_index_name = 'slopeTypeIndex'
        self.attrs[_slope_type_index_name] = SUMMAWorkflow._slope_type_index(
            slope_value=1,
            elements=self.hru,
            element_name=self.dims['hru'])

        # 3. `hruId` and `gruId` values
        # Note:
        #     If needed, the hard-coded names here can be transformed to be read
        #     from the input objects.
        if self.verbose:
            logger.info("Adding `hruId` and `gruId` attributes")
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
        if self.verbose:
            logger.info("Adding `hru2gruId` attributes")
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
        if self.verbose:
            logger.info("Calculating and adding `latitude` and `logitude` "
                "attributes")
        coords_defs = {
            'latitude': 'centroid_x',
            'longitude': 'centroid_y',
        }
        centroids = utils._calculate_centroids(self.hru)
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
        if self.verbose:
            logger.info("Calculating and adding `area` attributes")
        areas = utils._calculate_polygon_areas(self.hru, target_area_unit='m^2')
        area_unit = str(areas['area'].pint.units)
        self.attrs['HRUarea'] = xr.DataArray(
            data=areas['area'].pint.magnitude,
            coords={'hru': areas[hru_fid].values},
            attrs={'unit': area_unit})

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
        #             stream * 2).
        #    Wouter: uses 30m constant value
        #
        # 9. `downHRUindex` values
        #    Martyn: for non-contiguous HRUs, assign downHRUindex the
        #            riparian, contiguous one
        #    Darri: sets the downHRUindex based on elevation data
        #    Mohamed: all constant 0 values
        #    Wouter: Set the downHRUindex based on elevation data
        if self.verbose:
            logger.info("Adding geospatial layers' attributes")

        # 7. `geospatial` layers
        # 7.1 `tan_slope` layer
        if self.verbose:
            logger.info("Adding `tan_slope` attributes")
        _slope_name = 'tan_slope'
        slope = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_slope_name],
            _slope_name,
            'mean',
            'hru')
        # 7.2 `contourLength` layer
        if self.verbose:
            logger.info("Adding `contourLength` attributes")
        _contour_name = 'contourLength'
        contour = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_contour_name],
            _contour_name,
            'length',
            'hru')
        # 7.3 `downHRUindex` layer
        if self.verbose:
            logger.info("Adding `downHRUindex` attributes")
        _hruidx_name = 'downHRUindex'
        hru_index = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_hruidx_name],
            _hruidx_name,
            'index',
            'hru')
        # 7.4 `elevation` layer
        if self.verbose:
            logger.info("Adding `elevation` attributes")
        _elv_name = 'elevation'
        elv = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_elv_name],
            _elv_name,
            'mean',
            'hru')
        # 7.5 `vegTypeIndex` layer
        if self.verbose:
            logger.info("Adding `vegTypeIndex` attributes")
        _veg_name = 'vegTypeIndex'
        veg = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_veg_name],
            _veg_name,
            'majority',
            'hru')
        # 7.6 `soilTypeIndex` layer
        if self.verbose:
            logger.info("Adding `soilTypeIndex` attributes")
        _soil_name = 'soilTypeIndex'
        soil = SUMMAWorkflow._geolayer_info(
            self.geospatial_data[_soil_name],
            _soil_name,
            'majority',
            'hru')

        _geolayers = {
            _slope_name: slope,
            _contour_name: contour,
            _hruidx_name: hru_index,
            _elv_name: elv,
            _veg_name: veg,
            _soil_name: soil
        }

        self._geolayers = _geolayers

        # Fill invalid (<=0) values with `self.fillna` values in GeoLayers
        for key in self.fillna['geospatial_data'].keys():
            if key in _geolayers:
                _geolayers[key] = self._fill_na(
                    _geolayers[key],
                    self.fillna['geospatial_data'][key])

        # Adding Geospatial layers to `self.attrs`
        self.attrs.update(_geolayers)

        if self.verbose:
            logger.info("Adding local and global attributes of the Dataset")
        # Updating local attributes
        # Assign attributes to each variable
        for var_name, attrs in attributes_local_attrs_default.items():
            if var_name in self.attrs:
                self.attrs[var_name].attrs.update(attrs)
        # Updating global attributes
        self.attrs = self.attrs.assign_attrs(
            attributes_global_attrs_default
        )

        # Saving if instructed
        self._save_ds(
            ds=self.attrs,
            save=save,
            save_path=save_path
        )

        if self.verbose:
            logger.info("SUMMA attributes initialized successfully.")

        if return_ds:
            return self.attrs
        else:
            return

    def init_forcing(
        self,
        return_ds: bool = False,
        save: bool = False,
        save_nc_path: Optional[PathLike | str] = None,
        save_list_path: Optional[PathLike | str] = None,
    ) -> Optional[xr.Dataset]:
        """"Prepare forcing dataset for the SUMMA setup.

        The preparation involves:
            - Name changes and unit adjustments
            - Sorting (SUMMA-specific) NetCDF files as forcing inputs
            - Creating a text file listing all inputs

        Parameters
        ----------
        return_ds : bool, optional
            If True, returns the prepared dataset (default=False).

        Notes
        -----
        - Timezone naming scheme (TZ) follows IANA's convention found at the
          following link (version 2025b):
            https://data.iana.org/time-zones/releases/tzdb-2025b.tar.lz
        - Only NetCDF files are supported as inputs.
 
        """
        # Logging
        if self.verbose:
            logger.info("Initializing attributes for SUMMA workflow...")

        # If instructed to save, paths are necessary
        if save:
            if save_nc_path is None or save_list_path is None:
                raise ValueError("Missing paths to save files.")
        # Modify the init repr 
        self.auxillary['init'].append('forcing')

        # If the timezone is not provided, hard-coded assumption to `local`
        if self._forcing_attrs['forcing_time_zone']:
            _forcing_tz = self._forcing_attrs['forcing_time_zone'].lower()
        else:
            _forzing_tz = 'local'

        if self._forcing_attrs['target_time_zone']:
            _target_tz = self._forcing_attrs['target_time_zone'].lower()
        else:
            _target_tz = 'local'

        if self.verbose:
            logger.info("Assigning timezone")
        # Specify the `forcing` and `target` timezones and also the string
        # for the "fileManager" in the axuillary dictionary
        forcing_tz, target_tz, tz_info = SUMMAWorkflow._specify_tz(_forcing_tz, _target_tz)
        self.auxillary['tz_info'] = tz_info

        # File generation is done within this function as SUMMA does not have
        # other mechanisms in reading forcing files

        # Create the forcing dataset list file
        if save:
            list_file = pathlib.Path(save_list_path)

            # Create parent directories if they don't exist
            list_file.parent.mkdir(parents=True, exist_ok=True)
            list_file.touch()
            list_file.write_text('')

        # FIXME: Considering Dask/joblib/others for parallelization
        for forcing in self._forcing:
            # Extract filename
            filename = os.path.basename(forcing)

            # log forcing files being processed
            if self.verbose:
                logger.info(f"Processing forcing file: {filename}")

            # First read the forcing file, using Xarray
            ds = xr.open_dataset(forcing)

            # Change object's `time` encoding
            ds.time.encoding = {}
            ds.time.encoding = SUMMAWorkflow._specify_time_encodings(ds.time)

            # Specify the `dt_init` or time-step period in seconds
            if 'dt_init' not in self.auxillary.keys():
                self.auxillary['dt_init'] = utils._freq_seconds(pd.infer_freq(ds.time))

            # Rename variables
            ds = ds.rename_vars(self.forcing_vars)
            # Only select variables included in `self.forcing_vars`
            ds = ds[[*self.forcing_vars.values()]]

            # Change units
            ds = SUMMAWorkflow._unit_change(
                ds,
                self.forcing_units,
                self.forcing_to_units,
                self._ureg)

            # Assure the order of dimensions are similar to that of self.attrs
            if hasattr(self, 'attrs'):
                ds = ds.reindex(dims=self.attrs.dims)
            else:
                warnings.warn('Dimensions not sorted, run'
                    ' `init_attrs(...)` method first.')

            # Change timezone and assing tz-naive datetime64[ns] timestamps
            if target_tz not in ('local'):
                ds = ds.assign_coords({
                   'time': ds.time.to_index().tz_localize(forcing_tz).tz_convert(target_tz).tz_localize(None)
                })
            # Updating local attributes
            # Assign attributes to each variable
            for var_name, attrs in forcing_local_attrs_default.items():
                if var_name in ds:
                    ds[var_name].attrs.update(attrs)
            # If `local` in self._forcing_attrs, update local attributes
            if 'local' in self._forcing_attrs:
                for var_name, attrs in self._forcing_attrs['local'].items():
                    if var_name in ds:
                        ds[var_name].attrs.update(attrs)
                    else:
                        warnings.warn(f"Variable `{var_name}` not found in "
                            "forcing dataset; skipping local attributes "
                            "update.")

            # Updating global attributes
            ds = ds.assign_attrs(forcing_global_attrs_default)
            # If `global` in self._forcing_attrs, update global attributes
            if 'global' in self._forcing_attrs:
                ds.attrs.update(self._forcing_attrs['global'])

            # Specify dt_step from the data
            data_step = utils._freq_seconds(pd.infer_freq(ds.time))
            # Add `data_step`
            ds['data_step'] = float(data_step)

            # Find the dimension other than 'time' and rename it to 'subbasin'
            dims_to_rename = {dim: 'subbasin' for dim in ds.dims if dim != 'time'}
            ds = ds.rename_dims(dims_to_rename)
            # Also change the variable name of `dims_to_rename` if such a variable
            # exists in the dataset
            for old_dim, new_dim in dims_to_rename.items():
                if old_dim in ds.variables:
                    ds = ds.rename_vars({old_dim: new_dim})

            # Saving if instructed
            if save:
                # Save files
                self._save_ds(
                    ds=ds,
                    save=save,
                    save_path=os.path.join(save_nc_path, filename),
                    unlimited_dims=['time'],
                    nc_format='NETCDF4_CLASSIC',
                )

                # Create a list of forcing files
                with list_file.open('a', encoding='utf-8') as f:
                    f.write(f'{filename}\n')

            # Close the dataset
            ds.close()

        # Warn user if files are not being saved
        if not save:
            warnings.warn('Forcing files processed without being saved.')

        if self.verbose:
            logger.info("Forcing dataset processed/initialized successfully.")

        if return_ds:
            # Just a representation---really rough return
            return xr.open_mfdataset(self._forcing)
        else:
            return

    def init_cold_state(
        self,
        return_ds: bool = False,
        save: bool = False,
        save_path: Optional[PathLike | str] = None,
    ) -> Optional[xr.Dataset]:
        """Creating self.coldstate and optionally returning it

        Parameters
        ----------

        Returns
        -------

        Raises
        ------


        """
        # Logging
        if self.verbose:
            logger.info("Initializing cold state")

        # Modify init repr string
        if 'cold_state' not in self.auxillary['init']:
            self.auxillary['init'].append('cold_state')
        else:
            del self.cold_state

        # local variables for easier access
        _layer_keys = ('layer', 'layers')
        _state_keys = ('state', 'states')

        # Check if necessary keys are provided
        for key in self.cold_state_data.keys():
            if key.lower() in _layer_keys:
                layers = self.cold_state_data[key.lower()]
            elif key.lower() in _state_keys:
                states = self.cold_state_data[key.lower()]
            else:
                raise IndexError(f"Invalid key `{key}` in `cold_state` dictionary."
                    " Check documentations at summaflow.readthedocs.io for "
                    "valid options.")
        layers = self.cold_state_data['layers']
        states = self.cold_state_data['states']

        # Check if the minimum information is provided
        ## In `layers`
        _layers_needed_keys = {'nsoil', 'nsnow'}
        _layers_provided_keys = set(k.lower() for k in layers.keys())
        if _layers_needed_keys != _layers_provided_keys:
            raise IndexError("Missing key in `cold_state` layers. "
                " Check documentations at summaflow.readthedocs.io for "
                "valid options.")
        ## & `states`
        _states_needed_keys = {'mlayerdepth'}
        _states_provided_keys = set(k.lower() for k in states.keys())
        # Find the intersection of needed keys and what has been provided
        if not _states_needed_keys & _states_provided_keys:
            raise IndexError("Missing/invalid keys in `cold_state` state variables "
                "to define `mLayerDepth`. Check documentations at "
                "summaflow.readthedocs.io for valid options.")

        # If the sum of `nSoil` and `nSnow` is not greater than 1, raise
        # an exception
        if sum([layers['nSoil'], layers['nSnow']]) < 1:
            raise ValueError("At least 1 layer of `nSoil` is needed.")
        if int(layers['nSoil']) < 1:
            raise ValueError("At least 1 layer of `nSoil` is needed.")

        # Check the dtype of layers
        _vars = ['nSoil', 'nSnow']
        for _v in _vars:
            if not isinstance(layers[_v], (int, float)):
                raise TypeError(f"`{_v}` value must be a number.")

        if self.verbose:
            logger.info("Calculating `iLayerHeight` from `mLayerDepth`")
        # `iLayerHeight` need be calculated based on `mLayerDepth`
        ilayerheight = [float(0)] + list(itertools.accumulate(states['mLayerDepth']))

        if self.verbose:
            logger.info("Calculating dimensions for cold state")
        # Calculating dimensions based on the inputs
        # FIXME: `nSoil` and `nSnow` should not be hardcoded below
        ## `midToto`: mid layer indices of all soil plus snow layers
        if self.verbose:
            logger.info("Calculating `midToto`")
        mid_toto = np.arange(int(layers['nSoil']) + int(layers['nSnow']))
        ## `midSoil`: mid layer indices of all soil layers
        if self.verbose:
            logger.info("Calculating `midSoil`")
        mid_soil = np.arange(int(layers['nSoil']))
        ## `ifcToto`: interfaces between all layers in the combined soil
        ##            and snow profile (including top and bottom)
        if self.verbose:
            logger.info("Calculating `ifcToto`")
        ifc_toto = np.arange(int(layers['nSoil']) + int(layers['nSnow']) + 1)
        ## `scalarv`: scalar variables and parameters (degenerate dimension)
        ##            hard-coded to 1
        if self.verbose:
            logger.info("Assigning `scalarv` to 1")
        scalarv = np.arange(1)

        # Defining basic dimensions and variables for the file
        gru_fid = self.topology_attrs.get('gru_fid')
        hru_fid = self.topology_attrs.get('hru_fid')
        cold_state_coords = {
            'hru': self.hru[hru_fid],
            'midSoil': mid_soil,
            'midToto': mid_toto,
            'ifcToto': ifc_toto,
            'scalarv': scalarv,
        }
        # Length of HRUs provided
        _hru_len = len(self.hru[hru_fid])
        # Variables for the cold state object;
        # FIXME: `nSoil` and `nSnow` should be also processed like the state
        #        variable
        if self.verbose:
            logger.info("Adding variables: `hruId`, `dt_init`, `nSoil`, and `nSnow`")
        cold_state_variables = {
            'hruId': ('hru', self.hru[hru_fid]),
            'dt_init': (('scalarv', 'hru'), [[self.auxillary['dt_init']] * _hru_len]),
            'nSoil': (('scalarv', 'hru'), [[layers['nSoil']] * _hru_len]),
            'nSnow': (('scalarv', 'hru'), [[layers['nSnow']] * _hru_len]),
        }

        # Since `iLayerHeight` can be calculated from `mLayerDepth`,
        # if not provided, by the user, we add it to `states` dictionary.
        if 'iLayerHeight' not in states:
            states['iLayerHeight'] = ilayerheight

        # Additional `state` variables added by users
        for var in states.keys():
            if self.verbose:
                logger.info(f"Adding variables: `{var}`")
            # make an array out of the input variables
            value = np.asarray(states[var])

            # Specify dimensions
            dim_one = SUMMAWorkflow._cold_state_dim(var)
            dims = (dim_one, 'hru')
            dims_lengths = [len(cold_state_coords[dim]) for dim in dims]
            if self.verbose:
                logger.info(f"    dimensions: {dims} with lengths {dims_lengths}")

            # if shape is zero, meaning a 0-dimensional array, make it
            # one dimensional
            if value.ndim == 0:
                value = np.array([value])
            elif value.ndim == 1 or value.ndim == 2:
                value = np.asarray(value)
            else:
                raise ValueError("Currently 1- or 2-dimensional array-like objects "
                    "are supported.")

            # Checking if the length of given array matches any of its
            # dimensions first
            matched = [v for v in dims_lengths if v in value.shape]
            if matched:
                matched_prod = math.prod(matched)
            else:
                matched_prod = 1

            repeated_value = np.repeat(value, math.prod(dims_lengths) / matched_prod)
            repeated_value = repeated_value.reshape(dims_lengths)
            repeated_value = repeated_value.astype(float)

            cold_state_variables[var] = (dims, repeated_value)

        # Creating an empty xarray.Dataset for the file
        if self.verbose:
            logger.info(f"Adding variables to the `cold_state` Dataset")
        self.cold_state = xr.Dataset(
            data_vars=cold_state_variables,
            coords=cold_state_coords,
        )

        if self.verbose:
            logger.info(f"Updating `cold_state` Dataset's local and global attributes")
        # Updating local attributes
        # Assign attributes to each variable
        for var_name, attrs in cold_state_local_attrs_default.items():
            if var_name in self.cold_state:
                self.cold_state[var_name].attrs.update(attrs)
        # Updating global attributes
        self.cold_state = self.cold_state.assign_attrs(
            cold_state_global_attrs_default
        )

        # Assure the order of dimensions are similar to that of self.attrs
        if self.verbose:
            logger.info(f"Reordering `cold_state` hrus and grus to match dimensions of `attrs`")
        self.cold_state = self.cold_state.reindex(dims=self.attrs.dims)

        # Saving if instructed
        self._save_ds(
            ds=self.cold_state,
            save=save,
            save_path=save_path
        )

        if return_ds:
            return self.cold_state
        else:
            return

    def init_template(
        self,
        save: bool = False,
        save_path: Optional[PathLike | str] = None,
    ) -> None:
        """Copying template files for SUMMA runs"""
        if self.verbose:
            logger.info("Preparing template files to copy if instructed")
        if save:
            os.makedirs(save_path, exist_ok=True)
            if self.verbose:
                logger.info(f"Target path is: {save_path}")
 
        # Using the package root to find the data directory
        data_dir_path = self.package_root.joinpath('data')

        for item in data_dir_path.iterdir():
            if item.is_file():
                with as_file(item) as src_path:
                    if self.verbose:
                        logger.info(f"Copying file: {item.name}")
                    shutil.copy2(src_path, os.path.join(save_path, item.name))

        if self.verbose:
            logger.info("Template files initialized successfully.")

        return

    def init_trial(
        self,
        return_ds: bool = False,
        save: bool = False,
        save_path: Optional[PathLike | str] = None,
    ) -> Optional[xr.Dataset]:
        """Preparing trialParams.nc file for SUMMA setups.
        This should ideally facilitate calibration, so the progress made on
        this front has not been included here."""
        if self.verbose:
            logger.info("Preparing trial parameters Dataset")
        # Defining basic dimensions and variables for the file
        hru_fid = self.topology_attrs.get('hru_fid')
        gru_fid = self.topology_attrs.get('gru_fid')

        if self.verbose:
            logger.info("Adding `hru` and `gru` dimensions to the Dataset")
        # Coordinate variables
        trial_params_coords = {
            'hru': self.hru[hru_fid],
            'gru': self.gru[gru_fid],
        }
        if self.verbose:
            logger.info("Adding `hruId` and `gruId` variables to the Dataset")
        # Variables for the trial object;
        trial_params_variables = {
            'hruId': ('hru', self.hru[hru_fid]),
            'gruId': ('hru', self.gru[gru_fid]),
        }

        # Trial xarray.Dataset object
        self.trial_params = xr.Dataset(
            data_vars=trial_params_variables,
            coords=trial_params_coords,
        )

        if self.verbose:
            logger.info(f"Updating `trial_params` Dataset's local and global attributes")
        # Updating local attributes
        # Assign attributes to each variable
        for var_name, attrs in trial_params_local_attrs_default.items():
            if var_name in self.trial_params:
                self.trial_params[var_name].attrs.update(attrs)
        # Updating global attributes
        self.trial = self.trial_params.assign_attrs(
            trial_params_global_attrs_default
        )

        if self.verbose:
            logger.info(f"Reordering `trial_params` hrus and grus to match dimensions of `attrs`")
        # Assure the order of dimensions are similar to that of self.attrs
        self.trial = self.trial.reindex(dims=self.attrs.dims)

        # If saving instructed
        if save:
            self._save_ds(
                ds=self.trial_params,
                save=save,
                save_path=save_path
            )

        if self.verbose:
            logger.info("Trial parameters Dataset initialized successfully.")
        if return_ds:
            return self.trial_params
        else:
            return

    def init_decisions(
        self,
        return_dict: bool = False,
        save: bool = False,
        save_path: Optional[PathLike | str] = None,
    ) -> Optional[pd.DataFrame]:
        """Preparing modelDecisions.txt file for SUMMA setups."""
        if self.verbose:
            logger.info("Preparing model decisions")
        # Jinja2 template for model decisions
        template = self.jinja2_env.get_template("SUMMA_model_decisions_template.jinja")

        # Reading the default model decisions
        if self.verbose:
            logger.info("Populating default model decisions")

        # Temporary decision dictionary
        _decisions = model_decisions_default.copy()

        # If the decisions are provided, update the default ones
        if self.decisions is not None:
            for key, value in _decisions['models'].items():
                if key in self.decisions:
                    if self.decisions[key] != _decisions['models'][key]:
                        # Logging
                        if self.verbose:
                            logger.info(f"    Modifying decision `{key}` with `{self.decisions[key]}`")
                        # Exchange default decision with user-provided one
                        _decisions['models'][key] = self.decisions[key]

            for key in self.decisions.keys():
                if key not in _decisions['models']:
                    warnings.warn(f"Invalid key `{key}` in `decisions` dictionary. "
                                   "Check the documentation for valid keys.")
        else:
            # Logging
            if self.verbose:
                logger.info("Using default model decisions as no user-provided ones")

        # exchanging user-provided decisions with the populated _decisions
        # object
        self.decisions = _decisions

        # create content
        if self.verbose:
            logger.info("Rendering model decisions template")
        content = template.render(
            decisions_dict=self.decisions,
            version=__version__,
            models='models',
            comments='comments',
        )

        if save:
            pathlib.Path(save_path).write_text(content, encoding='utf-8') 

        if self.verbose:
            logger.info("Model decisions rendered successfully.")

        if return_dict:
            return self.decisions

    def run(
        self,
        save: bool = True,
        path: Optional[PathLike | str] = None,
        **kwargs: Optional[Any],
    ) -> None:
        """Preparing fileManager.txt file for SUMMA setups.

        Parameters
        ----------
        path: PathLike or str, optional
            Path to save the fileManager.txt file. If not provided, uses
            the model_path from settings.
        kwargs: dict, optional
            Additional key-value pairs to be included in the fileManager.txt.
            Valid keys include:
                - controlVersion
                - simStartTime
                - simEndTime
                - tmZoneInfo
                - outFilePrefix
                - settingsPath
                - forcingPath
                - outputPath
                - initConditionFile
                - attributeFile
                - trialParamFile
                - forcingListFile
                - decisionsFile
                - outputControlFile
                - globalHruParamFile
                - globalGruParamFile
                - vegTableFile
                - soilTableFile
                - generalTableFile
                - noahmpTableFile

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If `path` is not a string or PathLike object.
        ValueError
            If `settings['start_date']` or `settings['end_date']` is not provided
            in the settings dictionary.
        KeyError
            If an invalid key is provided in `kwargs`.
        """
        if self.verbose:
            logger.info(f"Running SUMMA workflow")
        # If path is not provided, use the model_path
        if not path:
            path = self.settings['model_path']
        else:
            # If path is a PathLike object, convert it to string
            if isinstance(path, PathLike):
                path = str(path)

        # Make sure the path exists
        os.makedirs(path, exist_ok=True)

        # If the path is not a string, raise an error
        if not isinstance(path, (str, PathLike)):
            raise TypeError("`path` argument must be a string or PathLike object.")

        # If invalid key is provided in `kwargs`, raise an error
        valid_keys = {
            'controlVersion', 'simStartTime', 'simEndTime', 'tmZoneInfo',
            'outFilePrefix', 'settingsPath', 'forcingPath', 'outputPath',
            'initConditionFile', 'attributeFile', 'trialParamFile',
            'forcingListFile', 'decisionsFile', 'outputControlFile',
            'globalHruParamFile', 'globalGruParamFile', 'vegTableFile',
            'soilTableFile', 'generalTableFile', 'noahmpTableFile'
        }
        if kwargs:
            invalid_keys = set(kwargs.keys()) - valid_keys
            if invalid_keys:
                raise KeyError(f"Invalid keys in `kwargs`: {invalid_keys}. "
                               "Check the documentation for valid keys.")

        # If `settings.['start_date']` or `settings['end_date']` or not provided,
        if 'start_date' not in self.settings or 'end_date' not in self.settings:
            raise ValueError("`settings['start_date']` and "
                             "`settings['end_date']` must be provided in the "
                             "settings dictionary.")

        # Initialize attributes
        self.auxillary['tz_info'] = SUMMAWorkflow._specify_tz(
            forcing=self._forcing_attrs['forcing_time_zone'],
            target=self._forcing_attrs['target_time_zone']
        )[-1]

        # Create the fileManager.txt content
        # key names are not PEP8 standard as they are used in SUMMA
        manager_dict = {
            'controlVersion': 'SUMMA_FILE_MANAGER_V3.0.0',
            'simStartTime': utils.format_date_string(self.settings['start_date']),
            'simEndTime': utils.format_date_string(self.settings['end_date']),
            'tmZoneInfo': self.auxillary['tz_info'],
            'outFilePrefix': 'run',
            'settingsPath': os.path.join(
                os.path.abspath(path),
                'settings',
                'SUMMA/'),
            'forcingPath': os.path.join(
                os.path.abspath(path),
                'forcing',
                'SUMMA/'),
            'outputPath': os.path.join(
                os.path.abspath(path),
                'output',
                'SUMMA/'),
            'initConditionFile': 'coldState.nc',
            'attributeFile': 'attributes.nc',
            'trialParamFile': 'trialParams.nc',
            'forcingListFile': 'forcingFileList.txt',
            'decisionsFile': 'modelDecisions.txt',
            'outputControlFile': 'outputControl.txt',
            'globalHruParamFile': 'localParamInfo.txt',
            'globalGruParamFile': 'basinParamInfo.txt',
            'vegTableFile': 'TBL_VEGPARM.TBL',
            'soilTableFile': 'TBL_SOILPARM.TBL',
            'generalTableFile': 'TBL_GENPARM.TBL',
            'noahmpTableFile': 'TBL_MPTABLE.TBL',
            }

        # If `kwargs` is provided, update the `manager_dict`
        if kwargs:
           for key, value in kwargs.items():
                if key in manager_dict:
                    manager_dict[key] = value
                else:
                    raise KeyError(f"Invalid key `{key}` in `kwargs`. "
                                   "Check the documentation for valid keys.")

        # Create the directories if they do not exist
        dirs = ('settingsPath', 'forcingPath', 'outputPath')
        for dir_key in dirs:
            dir_path = manager_dict[dir_key]
            if not os.path.exists(dir_path):
                os.makedirs(dir_path, exist_ok=True)

        # Write filemanager.txt content
        template = self.jinja2_env.get_template("SUMMA_file_manager_template.jinja")
        content = template.render(
            manager_dict=manager_dict,
            version=__version__,
        )

        if save:
            # 1. initialize attributes
            self.init_attrs(save=False)
            _ds = self.attrs.drop_indexes(
                ['hru', 'gru']  # drop indexes to avoid errors in SUMMA
            )
            _ds = _ds.reset_coords(
                ['hru', 'gru'],  # reset coordinates to avoid errors in SUMMA
                drop=True,
            )
            self._save_ds(
                ds=_ds,
                save=save,
                save_path=os.path.join(
                    manager_dict['settingsPath'],
                    'attributes.nc'),
            )

            # 2. initialize forcing
            self.init_forcing(
                save=save,
                save_nc_path=manager_dict['forcingPath'],
                save_list_path=os.path.join(
                    manager_dict['settingsPath'],
                    'forcingFileList.txt'),
            )

            # 3. initialize cold state
            self.init_cold_state(save=False)
            _ds = self.cold_state.drop_indexes(
                ['hru', 'midToto', 'midSoil', 'ifcToto', 'scalarv'],
            )
            _ds = _ds.reset_coords(
                ['hru', 'midToto', 'midSoil', 'ifcToto', 'scalarv'],
                drop=True,
            )
            self._save_ds(
                ds=_ds,
                save=save,
                save_path=os.path.join(
                    manager_dict['settingsPath'],
                    'coldState.nc'),
            )

            # 4. initialize trial parameters
            self.init_trial(save=False)
            _ds = self.trial_params.drop_indexes(
                ['hru', 'gru'],
            )
            _ds = _ds.reset_coords(
                ['hru', 'gru'],
                drop=True,
            )
            self._save_ds(
                ds=_ds,
                save=save,
                save_path=os.path.join(
                    manager_dict['settingsPath'],
                    'trialParams.nc'),
            )

            # 5. initialize model decisions
            self.init_decisions(
                save=save,
                save_path=os.path.join(
                    manager_dict['settingsPath'],
                    'modelDecisions.txt'),
            )
            # 6. initialize output control
            self.init_template(
                save=save,
                save_path=os.path.join(
                        manager_dict['settingsPath']),
            )
            # 7. initialize file manager
            pathlib.Path(
                os.path.join(
                    manager_dict['settingsPath'],
                    'fileManager.txt')
                ).write_text(content, encoding='utf-8') 

        return

    def _save_ds(
        self,
        ds,
        save,
        save_path,
        unlimited_dims: Sequence[str] = None,
        nc_format: str = 'NETCDF4',
    ) -> None:
        """Save the dataset to the defined path"""
        # Save the file
        if save:
            # Logging
            if self.verbose:
                logger.info(f"Saving dataset to {save_path}")

            if isinstance(save_path, (PathLike, str)):
                # Create the directory
                if utils.is_file(save_path, {'.nc', '.nc4'}):
                    os.makedirs(os.path.dirname(save_path), exist_ok=True)
                    ds.to_netcdf(
                        save_path,
                        format=nc_format,
                        unlimited_dims=unlimited_dims,
                    )
                else:
                    raise ValueError("File name missing.")

            else:
                raise TypeError("`save_path` argument must be a PathLike"
                    " or string object.")
        else: # for legibility
            pass

    def _fill_na(
        self,
        arr: np.array,
        fill_value: float,
    ) -> np.array:
        """Fill invalid values in an array with a specified fill value.

        Parameters
        ----------
        arr : np.arraylike
            The input array to fill invalid values.
        fill_value : float
            The value to replace invalid values with.

        Returns
        -------
        np.arraylike
            The array with invalid values replaced by the fill value.
        """
        arr[arr <= 0] = fill_value

        return arr

    # static methods
    @staticmethod
    def _cold_state_dim(
        var: str,
    ) -> str:
        """Specify the cold state file dimensions to be used for the selected
        state variable of the model

        Parameters
        ----------
        var: str
            The state variable name

        Returns
        -------
        dim: str
            The proper first dimension name of the state variable. The second
            dimension is always `hru`, therefore, has not been processed here.

        Raises
        ------
        TypeError
            If `var` is not of type ``str``
        """
        if not isinstance(var, str):
            raise TypeError(f'Invalid `{var}` variable type. It must be of '
                'data type string.')
        if var.lower() == 'mlayermatrichead':
            return 'midSoil'
        if var.lower().startswith('mlayer'):
            return 'midToto'
        elif var.lower().startswith('ilayer'):
            return 'ifcToto'
        elif var.lower().startswith('scalar') or var.lower() == 'nsoil' or (var.lower() == 'nsnow' or var.lower() == 'dt_init'):
            return 'scalarv'
        else:
            raise ValueError

    @staticmethod
    def _unit_change(
        ds: xr.Dataset,
        units: Dict[str, str],
        to_units: Dict[str, str] = None,
        missing_unit: str = 'dimensionless',
        unit_registry: pint.UnitRegistry = None,
    ) -> Tuple[xr.Dataset, pint.UnitRegistry]:
        """Changing units for an xarray.Dataset object

        Parameters
        ----------
        ds : |Dataset|
            A xarray.Dataset object containing at least one variable.
        units : dict of str keys and values
            The keys shows `ds` variables and corresponding values represent
            physical units associated with variables. If any variable is
            missing, it is assigned as 'dimensionless'.
        to_units : dict of str keys and values [defaults to `None`]
            Similar structure to `units` but representing target units of each
            variable.
        missing_unit : str [defaults to 'dimensionless']
            Unit value for variables defined in `ds` but not available in
            `units`. The behaviour defaults to 'dimensionless' unit.
        unit_registry : pint.UnitRegistry [default to ``None``]
            Pint unit registry to query physical units.

        Raises
        ------
        IndexError
            If |Dataset| does not have at least one variable.
        """
        # Check to see if the Dataset provides at least one variable
        if len(ds.variables) == 0:
            raise IndexError("`ds` must at least contain one variable.")

        # Check if all units provided are present in the forcing file(s)
        for k in units:
            if k not in ds:
                raise ValueError(f"item {k} defined in "
                                 "`forcing_unit_mapping` cannot be found"
                                 " in `forcing_name_mapping` values.")

        # If all elements of `variables` not found in `units`,
        # assign them to None
        for v in ds:
            if v not in units:
                units[v] = 'dimensionless'

        # Now assign the units
        ds = ds.pint.quantify(units=units, unit_registry=unit_registry)

        # If `to_units` is defined
        if to_units:
            ds = ds.pint.to(units=to_units)

        # Print the netCDF file
        ds = ds.pint.dequantify()

        return ds

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
            if forcing in ('utc', 'gmt') or forcing in time_zones:
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
        time_stamps: Sequence[np.datetime64],
    ) -> Dict[str, Dict[str, str]]:
        """Necessary adhoc modifications on the forcing object's
        encoding
        """
        # estimate the frequency offset value
        _freq = pd.infer_freq(time_stamps)
        # get the full name
        _freq_long = utils._freq_longname(_freq)

        # Encoding dictionary appropriate to be included as a local attribute;
        # The default starting date is hard-coded, as it is a safe date to
        # include;
        # Time string format is non-standard
        _encoding = {
            'units': f'{_freq_long} since 1900-01-01 00:00',
            'calendar': 'gregorian',
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

        if isinstance(layer.stats[stat_names], numbers.Real):
            # If stat_names is a single value, convert it to a list
            stat_values = [layer.stats[stat_names]]
        else:
            # If stat_names is a sequence, use it as is
            stat_values = layer.stats[stat_names]

        da = xr.DataArray(
            data=stat_values,
            coords={
                dim_name: layer.stats.data.index.values,
            },
            dims=[dim_name],
            attrs={'units': str(layer.unit)},
            name=layer_name,
        )

        return da

    @staticmethod
    def _is_valid_integer(s):
        try:
            int(s)
            return True
        except ValueError:
            return False
