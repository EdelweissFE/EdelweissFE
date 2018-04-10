# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:23:14 2015

@author: c8441141
"""
from distutils.core import setup, Extension
from Cython.Build import cythonize
from os.path import expanduser, join
import numpy 

directives = {'boundscheck':            False,
              'wraparound':             False,
              'nonecheck' :             False,
              'initializedcheck' :      False} # efficient access of memoryviews

# 1) rootDirectory where the libraries are located
rootDirectory = expanduser("~/constitutiveModelling/")

"""
Build Extension for the UEL base element, linked to the bftUserLibrary
"""
extensions = [Extension("*",
                sources = ["fe/elements/uelbaseelement/element.pyx"],
                        include_dirs=[join(rootDirectory,'bftUserLibrary', "include"), numpy.get_include()],
                         libraries= ['bftUserLibrary'],
                                 library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                         runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                        language='c++',)
                        ]  

"""
Build Extension for the UMAT material library, linked to the bftUserLibrary
"""

extensions += [Extension("*", ["fe/materials/umatlibrary.pyx"],                                    
                                 include_dirs=[numpy.get_include()] + [join(rootDirectory,'bftUserLibrary', "include")], 
                                 library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                                 runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                                 libraries= ['bftUserLibrary'],
                                 language="c++",)]
"""
Auxiliary extension modules
"""

extensions += [Extension("*", ["fe/utils/elementresultcollector.pyx"],                                    
                                 include_dirs=[numpy.get_include()],
                                 language="c++",)]

extensions += [Extension("*", ["fe/utils/csrgenerator.pyx"],                                    
                                 include_dirs=[numpy.get_include()],
                                 language="c++",)]


"""
Build The parallel NISTParallel solver with OpenMP
"""

extensions += [Extension("*",
                sources = ["fe/solvers/nonlinearimplicitstaticparallel.pyx"],
                        include_dirs=[join(rootDirectory,'bftUserLibrary', "include"), numpy.get_include()] ,
                        library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] , 
                        runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                         libraries= ['bftUserLibrary'],
                         language='c++',
                         extra_compile_args=['-fopenmp', '-Wno-maybe-uninitialized', ],
                         extra_link_args=['-fopenmp'],)
                        ]  

"""
Build The Pardiso Interface
"""

extensions += [Extension("*",
                sources = ["fe/linsolve/pardiso/pardiso.pyx", 'fe/linsolve/pardiso/pardisoInterface.cpp'],
                         include_dirs=[ numpy.get_include(), 'fe/linsolve/pardiso/include' ] ,
                         libraries= [ 'mkl_intel_thread', 'mkl_core', 'mkl_rt', 'iomp5', ],
                         language='c++',
                         )
                        ]  

"""
Compile!
"""
setup(ext_modules = cythonize(extensions,
                            compiler_directives=directives,
                            annotate=True
                            ))
