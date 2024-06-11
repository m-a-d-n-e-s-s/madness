# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "maddft"
copyright = "2024, Adrian Hurtado"
author = "Adrian Hurtado"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
    "sphinx.ext.mathjax",
]
latex_engine = 'xelatex'
latex_elements = {
    'preamble': r'''
    \usepackage{physics},
    '''
}
latex_additional_files = ['mra.sty']



myst_enable_extensions = [
    "dollarmath",  # Allows $ symbols for math expressions
    "amsmath",  # Enables support for amsmath LaTeX environments
    "deflist",  # Enables definition lists

    # add other extensions as needed
]


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "alabaster"
#html_theme = "pydata_sphinx_theme"
html_theme = "sphinx_rtd_theme"


#html_static_path = ["_static"]

bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "plain"
bibtex_reference_style = "label"
