"""Configuration file for the Sphinx documentation builder."""

from typing import Any

from docutils import nodes
from sphinx.application import Sphinx

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
    "sphinx_reredirects",
    "sphinx_tabs.tabs",
]

# Retired pages redirect to where their content lives. Targets are relative to the
# redirecting page.
redirects = {
    "documentation": "index.html",
    "events/index": "../index.html",
    "events/newsletter": "../overview.html#how-to-help",
    "faq": "overview.html#frequently-asked-questions",
    "guides/index": "../submissions.html",
}

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
    # The number of sections in the root toctree, so each stays in the bar instead of
    # collapsing into a "More" menu.
    "header_links_before_dropdown": 6,
    # Keeps the site bar and the theme's header stuck to the top of the window together.
    "sticky_banners": True,
    # One sidebar, on the right: the page's outline, then the pages in its section.
    "secondary_sidebar_items": ["page-toc", "section-nav"],
    "back_to_top_button": False,
    "primary_sidebar_end": [],
    "footer_start": ["ord-footer"],
    "footer_center": [],
    "footer_end": [],
}
# The viewer and the editor have no dark mode, so neither does the documentation.
html_context = {"default_mode": "light"}
html_show_sourcelink = False
# No left sidebar on wide screens; on narrow ones the theme keeps it as the drawer that
# holds the sections.
html_sidebars: dict[str, list[str]] = {"**": []}

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


def _section_root(app: Sphinx, pagename: str) -> str | None:
    """Returns the page under the root toctree that the given page sits beneath.

    Args:
        app: The Sphinx application, whose environment holds the toctrees.
        pagename: Name of the page, such as ``guides/templates``.

    Returns:
        The top-level page's name, which is ``pagename`` itself for a top-level page, or
        None for the root page and for a page no toctree includes.
    """
    parents = {}
    for parent, children in app.env.toctree_includes.items():
        for child in children:
            parents.setdefault(child, parent)
    page = pagename
    while parents.get(page) not in (None, app.config.root_doc):
        page = parents[page]
    if parents.get(page) != app.config.root_doc:
        return None
    return page


def _add_section_root(
    app: Sphinx,
    pagename: str,
    templatename: str,
    context: dict[str, Any],
    doctree: nodes.document | None,
) -> None:
    """Gives _templates/section-nav.html the title and link of the page's section.

    Args:
        app: The Sphinx application.
        pagename: Name of the page being rendered.
        templatename: The page's template; unused.
        context: The page's template context, which gains ``ord_section_root``.
        doctree: The page's doctree; unused.
    """
    root = _section_root(app, pagename)
    if root is None or root not in app.env.titles:
        return
    context["ord_section_root"] = {
        "title": app.env.titles[root].astext(),
        "url": context["pathto"](root),
        "current": root == pagename,
    }


def setup(app: Sphinx) -> None:
    """Registers the hooks above.

    Args:
        app: The Sphinx application.
    """
    # Runs before the theme's handler, which renders each sidebar template to drop empty
    # ones.
    app.connect("html-page-context", _add_section_root, priority=400)
