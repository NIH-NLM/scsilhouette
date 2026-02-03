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
    'sphinx.ext.autosummary',  # Add this
    'myst_parser',
]

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Mock imports for packages that might not be available during build
autodoc_mock_imports = []  # Add any problematic imports here if needed

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
