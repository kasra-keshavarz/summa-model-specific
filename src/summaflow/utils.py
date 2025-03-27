"""Common utilities needed for the SUMMA setup.

FIXME: Technically, the tools used here can be generalized for
       all models. Due to time limitation, some compromise were
       considered.
"""
# third-party libraries
import xarray as xr

# built-in imports
from typing import (
    Dict,
    Sequence,
)


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

