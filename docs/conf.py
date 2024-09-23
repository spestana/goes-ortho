# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

from goes_ortho import __version__

sys.path.insert(0, os.path.abspath("../src/goes_ortho"))

autosummary_generate = True

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "goes_ortho"
copyright = "2024, Steven Pestana"
author = "Steven Pestana"
version = __version__
release = __version__
# the page title will default to "<project> <relsease> documentation" e.g. goes_ortho 0.0.2 documentation

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "nbsphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# prevent execution of jupyter notebooks when building docs
nbsphinx_execute = "never"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = [
    "_static"
]  # this may generate a warning: WARNING: html_static_path entry '_static' does not exist, if we don't have anything in _static
