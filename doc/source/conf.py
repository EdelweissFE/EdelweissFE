# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# for execution python code in text

try:
    from StringIO import StringIO  # noqa: F401
except ImportError:
    pass

from importlib import import_module

from docutils import nodes
from docutils.parsers.rst import Directive  # noqa: F401
from pygments.lexer import RegexLexer
from pygments.token import Comment, Keyword, Literal, Name, Operator, Text
from sphinx import addnodes  # noqa: F401
from sphinx.directives.code import (  # noqa: F401
    CodeBlock,
    container_wrapper,
    dedent_lines,
)
from sphinx.highlighting import lexers

project = "EdelweissFE"
copyright = "2022, Matthias Neuner"
author = "Matthias Neuner"
release = "v22.07"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


# sys.path.insert(0, os.path.abspath("../../"))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
    "numpydoc",
]

templates_path = ["_templates"]
exclude_patterns = []

autosummary_generate = True
autoclass_content = "class"
autodoc_member_order = "groupwise"
# autodoc_typehints = "both"
# less crowded:
autodoc_typehints = "description"

autoclass_content = "init"

napoleon_use_admonition_for_notes = True
numpydoc_show_class_members = True
numpydoc_class_members_toctree = False
numpydoc_show_inherited_class_members = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_logo = "./logo.png"

html_static_path = ["_static"]

html_css_files = [
    "css/custom.css",
]


class EdelweissFELexer(RegexLexer):
    name = "EdelweissFE lexer"
    expression_in_quotes = r'(["\'])(?:(?=(\\?))\2.)*?\1'
    expression_no_quotes = r"\w+"
    expression_no_quotes = r'[^,\n\'"]+'
    equalsign_with_potential_whitespaces = r"\s*=\s*"
    tokens = {
        "root": [
            (r"\s*\*{2}.*\n", Comment.Singleline),
            (r",", Text),
            (r"\*{1}[^,\n]*", Keyword),
            (expression_in_quotes + r"\s*(?==)", Name.Variable),
            (r"(?<==)\s*" + expression_in_quotes, Literal.String),
            (expression_no_quotes + r"\s*(?==)", Name.Variable),
            (r"(?<==)\s*" + expression_no_quotes, Literal.Number),
            (r"=", Operator.Word),
            (r"[^=,\n]+", Text),
        ],
    }


lexers["edelweiss"] = EdelweissFELexer(startinline=True)
pygments_style = "nord"


class PrettyPrintDirective(CodeBlock):
    has_content = True
    optional_arguments = 1
    required_arguments = 1

    def run(self):
        module_path, member_name = self.arguments[0].rsplit(".", 1)
        member_data = getattr(import_module(module_path), member_name)

        table = nodes.table(cols=2)
        group = nodes.tgroup()
        head = nodes.thead()
        body = nodes.tbody()

        if "caption" in self.options:
            title = nodes.title(text=self.options["caption"])
            table += title

        table += group
        group += nodes.colspec(colwidth=6)
        group += nodes.colspec(colwidth=6)
        group += head
        group += body

        row = nodes.row()
        row += nodes.entry("", nodes.paragraph("", nodes.Text("Option")))
        row += nodes.entry("", nodes.paragraph("", nodes.Text("Description")))
        head += row

        for key, val in member_data.items():
            row = nodes.row()
            row += nodes.entry("", nodes.literal(text=key))
            row += nodes.entry("", nodes.paragraph("", nodes.Text(val)))
            body += row

        return [
            table,
        ]


def doi_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    # rendered = nodes.Text(text)
    uri = "http://dx.doi.org/" + text
    ref = nodes.reference(rawtext, text, refuri=uri)
    return [nodes.literal("", "", ref)], []


def setup(app):
    app.add_directive("pprint", PrettyPrintDirective)

    app.add_role("doi", doi_role)
