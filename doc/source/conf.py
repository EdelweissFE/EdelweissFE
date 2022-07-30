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
html_logo = "./logo.png"

html_static_path = ['_static']

html_css_files = [ 'css/custom.css', ]


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
from sphinx.directives.code import CodeBlock, dedent_lines, container_wrapper

class PrettyPrintDirective(CodeBlock):

    has_content = True
    optional_arguments = 1
    required_arguments = 1

    option_spec = CodeBlock.option_spec

    def get_codeblock_node(self, code, language):
        """this is copied from sphinx.directives.code.CodeBlock.run

        it has been changed to accept code and language as an arguments instead
        of reading from self

        """
        document = self.state.document

        literal = nodes.literal_block(code, code)
        literal['language'] = language

        caption = self.options.get('caption', '')
        if caption:
            try:
                literal = container_wrapper(self, literal, caption)
            except ValueError as exc:
                return [document.reporter.warning(text_type(exc), line=self.lineno)]

        self.add_name(literal)

        return [literal]

    def run(self):

        module_path, member_name = self.arguments[0].rsplit('.', 1)
        member_data = getattr(import_module(module_path), member_name)
        code = pformat(member_data, 2, width=80)

        cb = self.get_codeblock_node(code, 'python')

        return cb


def setup(app):
    app.add_directive('pprint', PrettyPrintDirective)
