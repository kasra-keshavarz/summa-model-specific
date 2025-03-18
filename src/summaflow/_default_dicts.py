"""
Default dictionaries
"""

from summaflow import __about__


attributes_global_attrs_default = {
    'author': 'University of Calgary',
    'license': 'GNU General Public License v3 (or any later version)',
    'purpose': (
        'Created by SUMMA automated configuration scripts (SUMMAFlow) ',
        'f{__about__.version}'
        ),
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
        },
    'latitude': {
        'long_name': 'Latitude of HRUs centroid',
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
