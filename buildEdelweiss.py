# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:23:14 2015

@author: c8441141
"""
from distutils.core import setup, Extension
from Cython.Build import cythonize
from os.path import expanduser, join, exists
import numpy 

directives = {'boundscheck': False,
              'wraparound':False,
              'nonecheck' : False}

gccFlags = ["-std=c++11", "-w", "-g3"] 
msvcFlags = []

compFlags = gccFlags

"""
Definition of libraries, library dirs and additional aux. libraries
Non-existing libs. are automatically removed!
"""

# 0) rootDirectory where the libraries are located
rootDirectory = expanduser("~/constitutiveModelling")

# 1) aux. libraries
auxLibraries = {'bftMechanics' : join(rootDirectory, 'bftMechanics/lib'),}

# 2) cython extension modules for elements
cythonElements = {  'uelCPE4' :         ['bftMechanics'],
                    'uelCPE4R':         ['bftMechanics'],
                    'uelCPS4':          ['bftMechanics'], 
                    'uelCPS4nonLocal':  ['bftMechanics'],
                                        
                    'uelCPE8RNonLocal': ['bftMechanics'], 
                    'uelCPS8R':         ['bftMechanics'], 
                    'uelCPS8RNonLocal': ['bftMechanics'], 
                    }

# 3) cython extension module for umat material library
umats = [
#         "Meschke",
         "ModLeon",
         "ModLeonNonLocal",
#         "ModLeonPlaneStress",
#         "ModLeonPlaneStressWrapped",
#         "ModLeonAdaptive",
#         "ModLeonAnalytical",
#         "ModLeonSemiExplicit",
#         "ModLeonSemiExplicitAdaptive",
#         "ShotLeon",
#         "SchaedlichSchweiger",
#         "Schuetz",
#         "linearElasticSolidificationCreep",
#         "Unteregger",
        ]
umatAuxLibs = ['bftMechanics',]

"""
Build Extension for the UMAT material library
"""

#filtering
umats = [umat for umat in umats if exists(join(rootDirectory, umat))]     
umatAuxLibs = [ libName for libName in umatAuxLibs if exists(auxLibraries[libName]) ]
extensions = [Extension("*", ["fe/materials/umatlibrary.pyx"],                                    
                                 include_dirs=[numpy.get_include()] + [join(rootDirectory,umat, "include") for umat in umats], 
                                 library_dirs= [join(rootDirectory,umat, "lib") for umat in umats] 
                                     + [auxLibraries[lib] for lib in umatAuxLibs],
                                 libraries= [umat for umat in umats]+[lib for lib in umatAuxLibs],    
                                 language="c++",
                                 extra_compile_args=compFlags)]

"""
Build Extensions for cython elements
"""

# filtering
cythonElements= {el : elExtraLibs for el, elExtraLibs in cythonElements.items() if exists(join(rootDirectory, el))}
for el, elementExtraLibs in cythonElements.items():
    
    # no filtering, all aux. libs must exist!
    extraLibDirs = [auxLibraries[d] for d in elementExtraLibs] 
    libs = [el]  + elementExtraLibs

    extensions.append( Extension("*",
                            sources=[join("fe/elements", el.lower(), "element.pyx")],
                            language="c++",
                            extra_compile_args=compFlags,
                            include_dirs=[
                                join(rootDirectory,el, "include"),                 
                                          numpy.get_include()],
                            library_dirs= [join(rootDirectory,el, "lib")] + extraLibDirs,
                            libraries= libs
                            ))
"""
... And all remaining .pyx files
"""

extensions += [Extension("fe/*/*",
                sources = ["fe/*/*.pyx"],
                        include_dirs=[ numpy.get_include()],
                        language='c++',
                        extra_compile_args=['-fopenmp'],
                        extra_link_args=['-fopenmp'],)
                        ]
      	
"""
Compile!
"""
setup(ext_modules = cythonize(extensions,
                            compiler_directives=directives))
