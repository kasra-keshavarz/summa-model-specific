"""
Default dictionaries
"""

from .__about__ import (
    __version__
)


attributes_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'purpose': ("Created by SUMMA automated configuration scripts (SUMMAFlow)"
        f" {__version__}"),
}

attributes_local_attrs_default = {
    'hruId': {
        'long_name': 'Index of hydrological response unit (HRU)',
        },
    'gruId': {
        'long_name': 'Index of grouped response unit (GRU)',
        },
    'hru2gruId': {
        'long_name': 'Index of GRU to which the HRU belongs',
        },
    'downHRUindex': {
        'long_name': 'Index of downslope HRU (0 = basin outlet)',
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
        },
    'contourLength': {
        'long_name': 'Contour length of HRU',
        },
    'slopeTypeIndex': {
        'long_name': 'Index defining slop',
        },
    'soilTypeIndex': {
        'long_name': 'Index defining soil type',
        },
    'vegTypeIndex': {
        'long_name': 'Index defining vegetation type',
        },
    'mHeight': {
        'long_name': 'Measurement height above bare ground',
        },
}


cold_state_local_attrs_default = {
    'hruId': {
        'long_name': 'Index of hydrological response unit (HRU)',
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
    'scalarCanopyliq': {
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
        'long_name': '',
        'units': '',
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
    'purpose': (
        'Create a cold state .nc file for initial SUMMA runs; ',
        'Created by SUMMA automated configuration scripts (SUMMAFlow) ',
        'f{__about__.version}'
        ),
}


trial_params_local_attrs_default = {
    'hruId': {
        'long_name': 'Index of hydrological response unit (HRU)',
        },
}


trial_params_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'purpose': (
        'Create a trial parameter .nc file for initial SUMMA runs; ',
        'Created by SUMMA automated configuration package (SUMMAFlow) ',
        'f{__about__.version}'
        ),
}


forcing_local_attrs_default = {
    'time': {
        'long_name': 'time',
        'standard_name': 'time',
        'calendar': 'gregorian',
        },
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
    'Conventions': 'CF-1.6',
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'purpose': (
        'Created by SUMMA automated configuration package (SUMMAFlow) ',
        'f{__about__.version}'
        ),
}


default_dims = {
    'hru': 'hru',
    'gru': 'gru',
}
