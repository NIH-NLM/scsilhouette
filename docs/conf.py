# Configuration file for Sphinx documentation

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

project = 'scsilhouette'
copyright = '2026, NIH-NLM'
author = 'NIH-NLM'
release = '1.0.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Napoleon settings for Google/NumPy docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
