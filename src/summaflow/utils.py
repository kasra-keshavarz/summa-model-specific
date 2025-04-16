"""Common utilities needed for the SUMMA setup.

FIXME: Technically, the tools used here can be generalized for
       all models. Due to time limitation, some compromise were
       considered.
"""
# third-party libraries
import xarray as xr
import geopandas as gpd
import pint_pandas
# specifics from third-party libraries 
from shapely.geometry import Polygon
from pyproj import CRS
from pint import UnitRegistry

# built-in imports
from collections import (
    Counter,
)
from typing import (
    Dict,
    Sequence,
    Any,
    List,
)
import warnings


def _init_empty_ds(
    variables: Dict[str, Sequence],
    dims: Dict[str, str],
) -> xr.Dataset:
    """Create an empty |Dataset| object with *hru* and *gru* dimensions

    Parameters
    ----------
    variables : Dict
        The keys must be *gru* and *hru*, with the values being the
        corresponding 
    dims : Dict
        The dimension names for the *gru* and *hru* `variables`

    Returns
    -------
    |Dataset| : a |Dataset| representing `hru` and `gru` values and dimensions

    Raises
    ------
    TypeError
        If either variables or dims is not a dictionary
    ValueError
        If required keys are missing
    """
    # routine checks
    # validate input types
    if not isinstance(variables, dict):
        raise TypeError("variables must be a dictionary")
    if not isinstance(dims, dict):
        raise TypeError("dims must be a dictionary")


    # validate input
    if not all(k in variables for k in ['gru', 'hru']):
        raise ValueError("variables must contain both 'gru' and 'hru' keys")
    if not all(k in dims for k in ['gru', 'hru']):
        raise ValueError("dims must contain both 'gru' and 'hru' keys")

    # create coordinate variables
    coords = {
        dims['gru']: variables['gru'],
        dims['hru']: variables['hru']
    }

    # Create empty dataset with these coordinates
    ds = xr.Dataset(coords=coords)

    # Return the goddamn xarray.Dataset
    return ds


def _calculate_centroids(
    gdf: gpd.GeoDataFrame,
    target_crs: str = 'EPSG:4326',
    default_crs: str = 'EPSG:4326',
    calculation_crs: str = None,
    geographic_fallback: str = 'EPSG:3857',
) -> gpd.GeoDataFrame:
    """
    Calculate centroids with flexible CRS handling and projection options.
    
    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Input GeoDataFrame with polygon geometries.
    target_crs : str or pyproj.CRS, optional
        CRS for output coordinates (default: 'EPSG:4326').
    default_crs : str or pyproj.CRS, optional
        CRS to assume if none is provided (default: 'EPSG:4326').
    calculation_crs : str or pyproj.CRS or None, optional
        Specific CRS to use for centroid calculations. If None and input is
        geographic, will use appropriate local CRS (default: None).
    geographic_fallback : str or pyproj.CRS, optional
        CRS to use if input is geographic and no calculation_crs provided
        (default: 'EPSG:3857').

    Returns
    -------
    |GeoDataFrame|
        Original GeoDataFrame with centroid columns added. The column labels
        are *centroid_x* and *centroid_y*
    """
    
    gdf = gdf.copy()
    
    # Handle missing CRS
    if gdf.crs is None:
        warnings.warn(f"No CRS provided - assuming {default_crs}", UserWarning)
        gdf = gdf.set_crs(default_crs)
    
    original_crs = gdf.crs
    is_geographic = original_crs.is_geographic
    
    # Calculate centroids with proper projection if needed
    if is_geographic:
        if calculation_crs is None:
            # Use fallback CRS for calculation
            calculation_crs = geographic_fallback
            warnings.warn(
                f"Using {calculation_crs} for centroid calculations on geographic data. "
                "For better accuracy, specify a local projection using calculation_crs parameter.",
                UserWarning
            )
        
        # Project to calculation CRS and calculate centroids
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            projected = gdf.to_crs(calculation_crs)
            centroids = projected.geometry.centroid.to_crs(original_crs)
    else:
        centroids = gdf.geometry.centroid
    
    # Convert to target CRS if different
    if str(original_crs) != str(target_crs):
        centroids = centroids.to_crs(target_crs)
    
    # Add results to DataFrame
    gdf['centroid'] = centroids
    gdf['centroid_x'] = centroids.x
    gdf['centroid_y'] = centroids.y
    
    return gdf


def _calculate_polygon_areas(
    gdf: gpd.GeoDataFrame,
    target_area_unit: str = 'm^2',
    equal_area_crs: str = 'ESRI:54009'
) -> gpd.GeoDataFrame:
    r"""Calculate polygon areas in a GeoDataFrame with proper unit handling.
 
    The area calculation follows these steps:
 
    1. For geographic CRS (e.g., EPSG:4326), transforms to equal-area projection
    2. For projected CRS, uses native units if detectable
    3. Converts to target units using Pint's dimensional analysis

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Input GeoDataFrame containing polygon geometries. If CRS is undefined,
        assumes EPSG:4326 (WGS84) with coordinates in degrees (\degree).
    target_area_unit : str, optional
        Output unit for area values. Common options:
        - 'm\ :sup:`2`\' for square meters (default)
        - 'km\ :sup:`2`\' for square kilometers
        - 'mi\ :sup:`2`\' for square miles
        - 'acre' for acres
    equal_area_crs : str, optional
        Equal-area projection to use when input CRS is geographic. The default
        'ESRI:54009' (World Mollweide) provides areas in m\ :sup:`2`\.
        Alternative options include:
        - 'EPSG:3410' (NSIDC Polar Stereographic North)
        - 'EPSG:6933' (WGS 84 / NSIDC EASE-Grid 2.0 Global)

    Returns
    -------
    geopandas.GeoDataFrame
        Copy of input with added 'area' column containing Pint quantities.
        The units will match ``target_area_unit`` when possible.

    Notes
    -----
    - Area calculations for geographic CRS use:
      \[
      A = \int \int \sqrt{g} \, d\phi \, d\lambda
      \]
      where \(g\) is the determinant of the metric tensor for the projection.
    - Unit conversions use Pint's dimensional analysis:
      \[
      1\ \text{m}^2 = 10^{-6}\ \text{km}^2 = 2.47105 \times 10^{-4}\ \text{acres}
      \]

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import Polygon
    >>> gdf = gpd.GeoDataFrame(
    ...     geometry=[Polygon([(0,0),(1,0),(1,1),(0,1)])],
    ...     crs="EPSG:4326"
    ... )
    >>> result = calculate_polygon_areas(gdf, target_area_unit='km²')
    """

    # Initialize Pint unit registry
    ureg = UnitRegistry()

    # Make a copy of the input GeoDataFrame to avoid modifying the original
    result_gdf = gdf.copy()

    # Check if CRS is defined, default to EPSG:4326 if not
    if result_gdf.crs is None:
        wrn_msg = """No CRS defined for GeoDataFrame.
        Assuming EPSG:4326 (WGS84) for area calculation."""
        warnings.warn(wrn_msg)
        result_gdf.crs = 'EPSG:4326'

    # Transform to equal area CRS if original is geographic
    try: 
        equal_area_gdf = result_gdf.to_crs(equal_area_crs)
    except pyproj.exceptions.CRSError:
        raise ValueError(f"Failed to transform to equal area CRS: {equal_area_crs}. "
                             "Please provide a valid equal area CRS.")

    # Try to determine the unit from the CRS
    try:
        crs_unit = CRS(equal_area_gdf.crs).axis_info[0].unit_name
        if crs_unit in ('metre', 'meter'):
            area_unit = ureg('m^2')
        elif crs_unit == 'US survey foot':
            area_unit = ureg('survey_foot^2')
        elif crs_unit == 'foot':
            area_unit = ureg('foot^2')
        else:
            warnings.warn(f"CRS has linear unit '{crs_unit}' which is not automatically "
                         "mapped to a Pint unit. Area values will be unitless.")
    except (AttributeError, IndexError):
        warnings.warn("Could not determine units from CRS. Area values will be unitless.")

    # Calculate areas in the equal area CRS
    areas = equal_area_gdf.geometry.area.astype(f'pint[{area_unit}]')

    # Convert to target unit if areas have units
    if hasattr(areas, 'pint'):
        try:
            target_unit = ureg(target_area_unit)
            converted_areas = areas.pint.to(target_unit)
        except:
            warnings.warn(f"Could not convert to target unit '{target_area_unit}'. "
                          "Keeping original units.")
            converted_areas = areas
    else:
        warnings.warn("Could not determine units from CRS. Area values will be unitless.")

    # Add area column to the result GeoDataFrame
    result_gdf['area'] = converted_areas

    return result_gdf


def unique_dict_values(
    d: Dict[Any, Any]
) -> Dict[Any, Any]:
    """Extract keys from a dictionary whose values appear exactly once.

    This function identifies all keys in a dictionary whose values have
    exactly one occurrence across all key-value pairs.

    Parameters
    ----------
    d : Dict[Any, Any]
        The input dictionary to analyze. Keys and values can be of any hashable type.

    Returns
    -------
    Dict[Any, Any]
        A list of keys whose values appear exactly once in the dictionary.
        The order is not guaranteed (depends on Python's hash implementation).

    Examples
    --------
    >>> unique_dict_values({'a': 1, 'b': 2, 'c': 1, 'd': 3})
    {'b': '2', 'd': '3'}

    >>> unique_dict_values({'x': 'apple', 'y': 'banana', 'z': 'apple'})
    {'y': 'banana'}

    Notes
    -----
    - The function uses collections.Counter for efficient counting
    - For Python 3.7+, you can use dict instead of OrderedDict as insertion
      order is preserved by default
    - All dictionary values must be hashable types
    """
    # Count occurrences of all values
    value_counts = Counter(d.values())
 
    # Get the set of unique values (appear exactly once)
    unique_values = {value for value, count in value_counts.items() if count == 1}
 
    # Return keys that map to these unique values
    return {key: value for key, value in d.items() if value in unique_values}


def nonunique_dict_values(
    d: Dict[Any, Any]
) -> Dict[Any, Dict]:
    """Extract keys from a dictionary whose values appear more than once.

    This function identifies all keys in a dictionary whose values have
    multiple occurrences across key-value pairs.

    Parameters
    ----------
    d : Dict[Any, Any]
        The input dictionary to analyze. Keys and values must be hashable.

    Returns
    -------
    Dict[Any, Any]
        A list of keys whose values appear more than once in the dictionary.
        The order is not guaranteed.

    Examples
    --------
    >>> nonunique_dict_values({'a': 1, 'b': 2, 'c': 1, 'd': 3})
    {1: ['c', 'a']}

    >>> nonunique_dict_values({'x': 'apple', 'y': 'banana', 'z': 'apple', 'w': 'banana'})
    {'banana': ['y', 'w'], 'apple': ['x', 'z']}

    Notes
    -----
    - If a value appears multiple times, all keys mapping to it are included.
    - Uses `collections.Counter` for efficient counting.
    - Values must be hashable.
    """
    value_counts = Counter(d.values())

    non_unique_values = {value for value, count in value_counts.items() 
        if count >= 2 and value is not None}

    # Create a dictionary referring to repeated values as refer corresponding
    # keys in a list - avoiding nested comprehensive dictionaries for
    # legibility
    dict_of_non_unique_values = {}
    for value in non_unique_values:
        dict_of_non_unique_values.update({value: [k for k, v in d.items() if v == value]})

    return dict_of_non_unique_values

def _freq_longname(
    freq_alias: 'str'
) -> str:
    """Returning fullname of a offset alias based on pandas conventions.

    Paramters
    ---------
    freq_alias : str
        Time offset alias which is usually a single character to represent
        time interval frequencies, such as 'H' for 'hours'

    Returns
    -------
    str
        fullname of the time offset
    """
    # Routine checks
    if not isinstance(freq_alias, str):
        raise TypeError(f"Frequency value of \'{freq_alias}\' is not"
                        "acceptable")

    # Check common time frequency aliases
    if freq_alias in ('H', 'h'):
        return 'hours'
    elif freq_alias in ('T', 'min'):
        return 'minutes'
    elif freq_alias in ('S'):
        return 'seconds'
    elif freq_alias in ('L', 'ms'):
        return 'milliseconds'
    else:
        raise ValueError(f"frequency value \'{freq_alias}\' is not"
                         "acceptable")

    return
