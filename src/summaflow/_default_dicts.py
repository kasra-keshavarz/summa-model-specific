"""
Default dictionaries mainly to define local and global attributes of
NetCDF files (xarray.Dataset objects).
"""

from .__about__ import __version__


attributes_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'workflow': ("Created by SUMMA automated configuration scripts (SUMMAFlow)"
        f" version {__version__}"),
}


attributes_local_attrs_default = {
    'hruId': {
        'long_name': 'Index of hydrological response unit (HRU)',
        'units': 'dimensionless',
        },
    'gruId': {
        'long_name': 'Index of grouped response unit (GRU)',
        'units': 'dimensionless',
        },
    'hru2gruId': {
        'long_name': 'Index of GRU to which the HRU belongs',
        'units': 'dimensionless',
        },
    'downHRUindex': {
        'long_name': 'Index of downslope HRU (0 = basin outlet)',
        'units': 'dimensionless',
        },
    'longitude': {
       'long_name': 'Longitude of HRUs centroid',
       'standard_name': 'longitude',
       'unit': 'degrees_east',
        },
    'latitude': {
        'long_name': 'Latitude of HRUs centroid',
        'standard_name': 'latitude',
        'unit': 'degress_north',
        },
    'elevation': {
        'long_name': 'Mean HRU elevation',
        },
    'HRUarea': {
        'long_name': 'Area of HRU',
        },
    'tan_slope': {
        'long_name': 'Average tangent slope of HRU',
        'units': 'dimensionless',
        },
    'contourLength': {
        'long_name': 'Contour length of HRU',
        },
    'slopeTypeIndex': {
        'long_name': 'Index defining slop',
        'units': 'dimensionless',
        },
    'soilTypeIndex': {
        'long_name': 'Index defining soil type',
        'units': 'dimensionless',
        },
    'vegTypeIndex': {
        'long_name': 'Index defining vegetation type',
        'units': 'dimensionless',
        },
    'mHeight': {
        'long_name': 'Measurement height above bare ground',
        },
}


cold_state_local_attrs_default = {
    'hruId': {
        'long_name': 'Index of hydrological response unit (HRU)',
        'units': 'dimensionless',
        },
    'dt_init': {
        'long_name': 'Time step size of forcing data',
        'units': 'second',
        },
    'nSoil': {
        'long_name': 'Number of soil layers',
        'units': 'dimensionless',
        },
    'nSnow': {
        'long_name': 'Number of snow layers',
        'units': 'dimensionless',
        },
    'scalarCanopyIce': {
        'long_name': 'Mass of ice on the vegetation canopy',
        'units': 'kilogram / meter ** 2',
        },
    'scalarCanopyLiq': {
        'long_name': 'Mass of liquid water on the vegetation canopy',
        'units': 'kilogram / meter ** 2',
        },
    'scalarCanairTemp': {
        'long_name': 'Temperature of the canopy air space',
        'units': 'pascal',
        },
    'scalarCanopyTemp': {
        'long_name': 'Temperature of the vegetation canopy',
        'units': 'kelvin',
        },
    'scalarSnowAlbedo': {
        'long_name': 'Snow albedo for the entire spectral band',
        'units': 'dimensionless',
        },
    'scalarSnowDepth': {
        'long_name': 'Total snow depth',
        'units': 'meter',
        },
    'scalarSWE': {
        'long_name': 'Snow water equivalent',
        'units': 'kilogram / meter ** 2',
        },
    'scalarSfcMeltPond': {
        'long_name': 'Ponded water caused by melt of the "snow without a layer"	',
        'units': 'kilogram / meter ** 2',
        },
    'scalarAquiferStorage': {
        'long_name': 'Relative aquifer storage; above bottom of the soil profile',
        'units': 'meter',
        },
    'mLayerDepth': {
        'long_name': 'Depth of each layer',
        'units': 'meter',
        },
    'iLayerHeight': {
        'long_name': 'Height of the layer interface; top of soil is 0.',
        'units': 'meter',
        },
    'mLayerTemp': {
        'long_name': 'Temperature of each layer',
        'units': 'kelvin',
        },
    'mLayerVolFracIce': {
        'long_name': 'Volumetric fraction of ice in each layer',
        'units': 'dimensionless',
        },
    'mLayerVolFracLiq': {
        'long_name': 'Volumetric fraction of liquid water in each layer',
        'units': 'dimensionless',
        },
    'mLayerMatricHead': {
        'long_name': 'Matric head of water in the soil',
        'units': 'dimensionless',
        },
}


cold_state_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'workflow': (
        'Create a cold state .nc file for initial SUMMA runs; '
        'Created by SUMMA automated configuration scripts (SUMMAFlow) '
        f'version {__version__}'
        ),
}


trial_params_local_attrs_default = {
    'hruId': {
        'long_name': 'Index of hydrological response unit (HRU)',
        'units': 'dimensionless',
        },
    'gruId': {
        'long_name': 'Index of grouped response unit (GRU)',
        'units': 'dimensionless',
        },
}


trial_params_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'workflow': (
        'Create a trial parameter .nc file for initial SUMMA runs; '
        'Created by SUMMA automated configuration package (SUMMAFlow) '
        f'version {__version__}'
        ),
}


forcing_local_attrs_default = {
    'latitude': {
        },
    'longitude': {
        },
    'hruId': {
        },
    'airpres': {
        },
    'LWRadAtm': {
        },
    'SWRadAtm': {
        },
    'pptrate': {
        },
    'airtemp': {
        },
    'spechum': {
        },
    'windspd': {
        },
    'data_step': {
        },
}


forcing_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'workflow': (
        'Created by SUMMA automated configuration package (SUMMAFlow) '
        f'version {__version__}'
        ),
}


model_decisions_default = {
    'comments': {
        'soilCatTbl': '(03) soil-category dateset (STAS, STAS-RUC, ROSETTA)',
        'vegeParTbl': '(04) vegetation category dataset',
        'soilStress': '(05) choice of function for the soil moisture control on stomatal resistance',
        'stomResist': '(06) choice of function for stomatal resistance',
        'num_method': '(07) choice of numerical method',
        'fDerivMeth': '(08) method used to calculate flux derivatives',
        'LAI_method': '(09) method used to determine LAI and SAI',
        'f_Richards': '(10) form of Richard\'s equation',
        'groundwatr': '(11) choice of groundwater parameterization',
        'hc_profile': '(12) choice of hydraulic conductivity profile',
        'bcUpprTdyn': '(13) type of upper boundary condition for thermodynamics',
        'bcLowrTdyn': '(14) type of lower boundary condition for thermodynamics',
        'bcUpprSoiH': '(15) type of upper boundary condition for soil hydrology',
        'bcLowrSoiH': '(16) type of lower boundary condition for soil hydrology',
        'veg_traits': '(17) choice of parameterization for vegetation roughness length and displacement height',
        'canopyEmis': '(18) choice of parameterization for canopy emissivity',
        'snowIncept': '(19) choice of parameterization for snow interception',
        'windPrfile': '(20) choice of wind profile through the canopy',
        'astability': '(21) choice of stability function',
        'canopySrad': '(22) choice of canopy shortwave radiation method',
        'alb_method': '(23) choice of albedo representation',
        'compaction': '(24) choice of compaction routine',
        'snowLayers': '(25) choice of method to combine and sub-divide snow layers',
        'thCondSnow': '(26) choice of thermal conductivity representation for snow',
        'thCondSoil': '(27) choice of thermal conductivity representation for soil',
        'spatial_gw': '(28) choice of method for the spatial representation of groundwater',
        'subRouting': '(29) choice of method for sub-grid routing',
    },
    'models': {
        'soilCatTbl': 'ROSETTA',
        'vegeParTbl': 'MODIFIED_IGBP_MODIS_NOAH',
        'soilStress': 'NoahType',
        'stomResist': 'BallBerry',
        'num_method': 'itertive',
        'fDerivMeth': 'analytic',
        'LAI_method': 'monTable',
        'f_Richards': 'mixdform',
        'groundwatr': 'bigBuckt',
        'hc_profile': 'constant',
        'bcUpprTdyn': 'nrg_flux',
        'bcLowrTdyn': 'zeroFlux',
        'bcUpprSoiH': 'liq_flux',
        'bcLowrSoiH': 'drainage',
        'veg_traits': 'Raupach_BLM1994',
        'canopyEmis': 'difTrans',
        'snowIncept': 'lightSnow',
        'windPrfile': 'logBelowCanopy',
        'astability': 'louisinv',
        'canopySrad': 'BeersLaw',
        'alb_method': 'conDecay',
        'compaction': 'anderson',
        'snowLayers': 'CLM_2010',
        'thCondSnow': 'jrdn1991',
        'thCondSoil': 'funcSoilWet',
        'spatial_gw': 'localColumn',
        'subRouting': 'qInstant',
    },
}


default_dims = {
    'hru': 'hru',
    'gru': 'gru',
}
