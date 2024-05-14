Input file syntax
=================

Relevant module: ``edelweissfe.utils.inputfileparser.py``

EdelweissFE uses an input file syntax similar to Abaqus, using human readable input files.
The general syntax makes use of ``*keywords``, possibly with ``option=values`` in the same line and
supplemented by data lines.
In general, options and data entries are comma separated.
Indentations and blank lines have no effect.

.. code-block:: edelweiss

    *keyword, option1=value, option2=another_value
    data, data, data, data
    data, data, data, data
    data, data, data, data

    **this is a comment
    **again a comment

    *another_keyword, option1=value, option2=another_value
    data, data, data, data


The available input syntax is defined in ``fe/utils/inputfileparser.py``,
and available keywords can be printed using

.. code-block:: edelweiss

    python edelweiss.py --keywords

.. automodule:: edelweissfe.utils.inputfileparser
   :members:

Keywords
********

.. include:: keywords.rst

