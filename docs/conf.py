# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'summaflow'
copyright = '2022-2025, University of Calgary'
author = 'Kasra Keshavarz, Wouter Knoben, Ignacio Aguirre'
release = '0.0.1-dev0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

# -- Extra options -----------------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'geopandas': ('https://geopandas.org/en/stable/', None),
}

# -- Shortcuts ---------------------------------------------------------------
rst_prolog = """
.. |GeoDataFrame| replace:: :class:`geopandas.GeoDataFrame`
.. |PathLike| replace:: :class:`os.PathLike`
"""
