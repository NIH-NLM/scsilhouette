# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# Added after running the Sphinx configuration command:
# sphinx-quickstart docs/ --quiet -p "scsilhouette" -a "NIH/NLM" --sep --makefile --no-batchfile
#
# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'scsilhouette'
copyright = '2025, NIH/NLM'
author = 'NIH/NLM'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',   # for Google-style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'myst_parser',           # <-- for Markdown support
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

autosummary_generate = True  # auto-generate .rst files for modules


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
templates_path = ['_templates']

html_theme_options = {
    "navigation_with_keys": True,
    "navigation_depth": 4,
    "collapse_navigation": False,
    "canonical_url": "https://nih-nlm.github.io/scsilhouette.github.io/"
}

html_baseurl = "https://nih-nlm.github.io/scsilhouette.github.io"
html_extra_path = ['_static']
