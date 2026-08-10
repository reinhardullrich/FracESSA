from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

project = "FracESSA"
author = "Reinhard Ullrich"

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "breathe"]
autodoc_typehints = "description"
napoleon_google_docstring = True
breathe_projects = {"FracESSA": str(ROOT / "docs" / "_build" / "xml")}
breathe_default_project = "FracESSA"

html_theme = "alabaster"
html_title = "FracESSA Documentation"
html_baseurl = "https://reinhardullrich.github.io/fracessa/docs/"
html_theme_options = {
    "description": "Fractional ESS Analyzer",
    "github_user": "reinhardullrich",
    "github_repo": "fracessa",
    "github_button": True,
}
