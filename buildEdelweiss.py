# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:23:14 2015

@author: c8441141
"""
from distutils.core import setup, Extension
from Cython.Build import cythonize
from os.path import expanduser, join
import numpy 

directives = {'boundscheck':    False,
              'wraparound':     False,
              'nonecheck' :     False}

gccFlags = ["-std=c++11", "-w", "-g3"] 
msvcFlags = []

compFlags = gccFlags

"""
Definition of libraries, library dirs and additional aux. libraries
"""

# 1) rootDirectory where the libraries are located
rootDirectory = expanduser("~/constitutiveModelling")

# 2) cython extension modules for elements
bftUserLibraryElements = [  
                             'uelCPE4',
                             'uelCPE4R',
                             'uelCPS4',
                             'uelCPS4NonLocal',
                             'uelCPE8RNonLocal',
                             'uelCPS8R',
                             'uelCPS8RNonLocal',
                             'uelCPS8NonLocal',
                            ]

"""
Build Extension for the UMAT material library, linked to the bftUserLibrary
"""

extensions = [Extension("*", ["fe/materials/umatlibrary.pyx"],                                    
                                 include_dirs=[numpy.get_include()] + [join(rootDirectory,'bftUserLibrary', "include")], 
                                 library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                                 runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                                 libraries= ['bftUserLibrary'],
                                 language="c++",
                                 extra_compile_args=compFlags)]

"""
Build Extensions for UEL Libraries, linked to the bftUserLibrary
"""

for el in bftUserLibraryElements:
    
    extensions.append( Extension("*",
                        sources=[join("fe/elements", el.lower(), "element.pyx")] ,
                        language="c++",
                        extra_compile_args=compFlags,
                        include_dirs=[
                            join(rootDirectory,'bftUserLibrary', "include"), numpy.get_include()] ,
                         runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                         libraries= ['bftUserLibrary'],
                        ))
"""
The parallel NISTParallel solver with OpenMP
"""

extensions += [Extension("*",
                sources = ["fe/solvers/nonlinearimplicitstaticparallel.pyx"],
                        include_dirs=[join(rootDirectory,'bftUserLibrary', "include"), numpy.get_include()] ,
                         runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                         libraries= ['bftUserLibrary'],
                        language='c++',
                        extra_compile_args=['-fopenmp'],
                        extra_link_args=['-fopenmp'],)
                        ]  
    

extensions += [Extension("*",
                sources = ["fe/elements/AbstractBaseElements/CPPBackendedElement.pyx"],
                        include_dirs=[ numpy.get_include()]+['fe/utils/cythonElementBackends/include'],
                        language='c++',)
                        ]  
    
"""
... And all remaining .pyx files which require no special setup
"""
# extensions += [Extension("fe/*/*",
                # sources = ["fe/*/*.pyx"],
                        # include_dirs=[numpy.get_include()],
                        # language='c++',)]
      	
"""
Compile!
"""
setup(ext_modules = cythonize(extensions,
                            compiler_directives=directives))
