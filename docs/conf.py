from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

project = "FracESSA"
author = "Reinhard Ullrich"
copyright = "2026, Reinhard Ullrich"

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "breathe"]
autodoc_typehints = "description"
napoleon_google_docstring = True
breathe_projects = {"FracESSA": str(ROOT / "docs" / "_build" / "xml")}
breathe_default_project = "FracESSA"

html_theme = "alabaster"
html_title = "FracESSA"
html_short_title = "FracESSA"
html_baseurl = "https://reinhardullrich.github.io/fracessa/"
html_static_path = ["_static"]
html_css_files = ["fracessa.css"]
html_show_sphinx = False
html_sidebars = {"**": ["about.html", "searchfield.html", "navigation.html"]}
html_theme_options = {
    "description": "Exact ESS analysis for symmetric rational games",
    "fixed_sidebar": True,
    "github_user": "reinhardullrich",
    "github_repo": "fracessa",
    "github_button": False,
    "sidebar_collapse": False,
}
