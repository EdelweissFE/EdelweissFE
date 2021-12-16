# -*- coding: utf-8 -*-
#  ---------------------------------------------------------------------
#
#  _____    _      _              _         _____ _____
# | ____|__| | ___| |_      _____(_)___ ___|  ___| ____|
# |  _| / _` |/ _ \ \ \ /\ / / _ \ / __/ __| |_  |  _|
# | |__| (_| |  __/ |\ V  V /  __/ \__ \__ \  _| | |___
# |_____\__,_|\___|_| \_/\_/ \___|_|___/___/_|   |_____|
#
#
#  Unit of Strength of Materials and Structural Analysis
#  University of Innsbruck,
#  2017 - today
#
#  Matthias Neuner matthias.neuner@uibk.ac.at
#
#  This file is part of EdelweissFE.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  The full text of the license can be found in the file LICENSE.md at
#  the top level directory of EdelweissFE.
#  ---------------------------------------------------------------------
"""
Created on Thu May 21 14:23:14 2015

@author: c8441141
"""
from distutils.core import setup, Extension
from Cython.Build import cythonize
from os.path import expanduser, join
import numpy
import sys
import os

directives = {
    "boundscheck": False,
    "wraparound": False,
    "nonecheck": False,
    "initializedcheck": False,
}

default_install_prefix = sys.prefix
print("*" * 80)
print("EdelweissFE setup")
print("System prefix: " + sys.prefix)
print("*" * 80)

marmot_dir = expanduser(os.environ.get("MARMOT_INSTALL_DIR", default_install_prefix))
mkl_include = expanduser(os.environ.get("MKL_INCLUDE_DIR", join(default_install_prefix, "include")))
eigen_include = expanduser(os.environ.get("EIGEN_INCLUDE_DIR", join(default_install_prefix, "include")))
print("Marmot install directory (overwrite via environment var. MARMOT_INSTALL_DIR):")
print(marmot_dir)
print("MKL include directory (overwrite via environment var. MKL_INCLUDE_DIR):")
print(mkl_include)
print("Eigen directory (overwrite via environment var. EIGEN_INCLUDE_DIR):")
print(eigen_include)
print("*" * 80)

"""
Build Extension for the MarmotElement base element, linked to the Marmot library
"""
extensions = [
    Extension(
        "*",
        sources=["fe/elements/marmotelement/element.pyx"],
        include_dirs=[join(marmot_dir, "include"), numpy.get_include()],
        libraries=["Marmot"],
        library_dirs=[join(marmot_dir, "lib")],
        runtime_library_dirs=[join(marmot_dir, "lib")],
        language="c++",
    )
]
"""
Auxiliary extension modules
"""
extensions += [
    Extension(
        "*",
        ["fe/utils/elementresultcollector.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++",
    )
]

extensions += [
    Extension(
        "*",
        ["fe/utils/csrgenerator.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++",
    )
]

"""
Build The parallel NISTParallelMK2 solver with OpenMP
"""
extensions += [
    Extension(
        "*",
        sources=["fe/solvers/nonlinearimplicitstaticparallelmk2.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++",
        extra_compile_args=[
            "-fopenmp",
            "-Wno-maybe-uninitialized",
        ],
        extra_link_args=["-fopenmp"],
    )
]

"""
Build The (old) parallel NISTParallel solver with OpenMP
"""
extensions += [
    Extension(
        "*",
        sources=["fe/solvers/nonlinearimplicitstaticparallel.pyx"],
        include_dirs=[numpy.get_include()] + [join(marmot_dir, "include")],
        language="c++",
        extra_compile_args=[
            "-fopenmp",
            "-Wno-maybe-uninitialized",
        ],
        extra_link_args=["-fopenmp"],
    )
]

"""
Build The Pardiso Interface
"""
extensions += [
    Extension(
        "*",
        sources=[
            "fe/linsolve/pardiso/pardiso.pyx",
            "fe/linsolve/pardiso/pardisoInterface.cpp",
        ],
        include_dirs=[
            numpy.get_include(),
            mkl_include,
            eigen_include,
        ],
        libraries=[
            "mkl_intel_thread",
            "mkl_core",
            "mkl_rt",
            "iomp5",
        ],
        language="c++",
    )
]

"""
Compile!
"""
setup(ext_modules=cythonize(extensions, compiler_directives=directives, annotate=True, language_level=3))
