# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'EdelweissFE'
copyright = '2022, Matthias Neuner'
author = 'Matthias Neuner'
release = 'v22.07'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

extensions = [ 'sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = []
html_logo = "./logo.png"


# for execution python code in text
import sys
from os.path import basename

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from importlib import import_module
from pprint import pformat
from docutils.parsers.rst import Directive
from docutils import nodes
from sphinx import addnodes

class PrettyPrintDirective(Directive):
    """Render a constant using pprint.pformat and insert into the document"""
    required_arguments = 1

    def run(self):
        module_path, member_name = self.arguments[0].rsplit('.', 1)

        member_data = getattr(import_module(module_path), member_name)
        code = pformat(member_data, 2, width=80)
       
        literal = nodes.literal_block(code, code)
        literal['language'] = 'python'

        return [
                addnodes.desc_name(text=''),
                addnodes.desc_content('', literal)
        ]


def setup(app):
    app.add_directive('pprint', PrettyPrintDirective)
