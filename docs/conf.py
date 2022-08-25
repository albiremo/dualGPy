# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
from recommonmark.parser import CommonMarkParser

sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DualGPy'
copyright = '2022, Alberto Remigi'
author = 'Alberto Remigi'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc",
    "sphinx.ext.imgmath",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx_rtd_theme",
    "recommonmark",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = ['.rst', '.md']
source_parsers = {
    '.md': CommonMarkParser,
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
