"""Configuration file for the Sphinx documentation builder."""

# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute.

# -- Project information -----------------------------------------------------

project = "Open Reaction Database"
# (sphinx requires this exact name)
copyright = "2020–2026 Open Reaction Database Project Authors"  # noqa: A001
author = "Open Reaction Database Project Authors"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",
    "sphinx_tabs.tabs",
]

# A docstring's Attributes section renders as fields of its class, since autodoc already
# gives each attribute an entry of its own and two entries for one name collide.
napoleon_use_ivar = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The site bar above the theme's header (_templates/layout.html) carries the logo and
# the links shared with the viewer, so the theme's own header holds only the
# documentation's sections, one per entry in the root toctree, and the search field.
html_theme = "pydata_sphinx_theme"
html_title = "Open Reaction Database"
html_logo = "images/ord-logo.png"
html_theme_options = {
    "navbar_start": [],
    "navbar_center": ["navbar-nav"],
    "navbar_end": [],
    "navbar_persistent": ["search-button-field"],
    "navbar_align": "left",
    # Keeps the site bar and the theme's header stuck to the top of the window together.
    "sticky_banners": True,
    "secondary_sidebar_items": ["page-toc"],
    "back_to_top_button": False,
    "primary_sidebar_end": [],
    "footer_start": ["ord-footer"],
    "footer_center": [],
    "footer_end": [],
}
# The viewer and the editor have no dark mode, so neither does the documentation.
html_context = {"default_mode": "light"}
html_show_sourcelink = False
html_sidebars = {"**": ["sidebar-nav-bs"], "index": []}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    "https://fonts.googleapis.com/css2?family=Roboto:ital,wght@0,100..900;1,100..900&display=swap",
    "css/custom.css",
]

master_doc = "index"
