# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:23:14 2015

@author: c8441141
"""
from distutils.core import setup, Extension
from Cython.Build import cythonize
from os.path import expanduser, join, exists
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

# 0) rootDirectory where the libraries are located
rootDirectory = expanduser("~/constitutiveModelling")

# 1) aux. libraries
auxLibraries = {'bftMechanics' : join(rootDirectory, 'bftMechanics/lib'),}

# 2) cython extension modules for elements
cythonElements = {  
#                    'uelCPE4' :         ['bftMechanics'],
#                    'uelCPE4R':         ['bftMechanics'],
                    'uelCPS4':          ['bftUserLibrary'], 
#                    'uelCPS4NonLocal':  ['bftMechanics'],
                                        
#                    'uelCPE8RNonLocal': ['bftMechanics'], 
#                    'uelCPS8R':         ['bftMechanics'], 
#                    'uelCPS8RNonLocal': ['bftMechanics'], 
#                    'uelCPS8NonLocal':  ['bftMechanics'], 
                    }

# 3) cython extension module for umat material library
#umats = [
##        "linearElastic",
##         "Meschke",
#         "ModLeon",
#         "ModLeonNonLocal",
##         "ModLeonPlaneStress",
##         "ModLeonPlaneStressWrapped",
##         "ModLeonAdaptive",
##         "ModLeonAnalytical",
##         "ModLeonSemiExplicit",
##         "ModLeonSemiExplicitAdaptive",
##         "ShotLeon",
##         "SchaedlichSchweiger",
##         "Schuetz",
##         "linearElasticSolidificationCreep",
##         "Unteregger",
#        ]
#umatAuxLibs = ['bftMechanics',]

"""
Build Extension for the UMAT material library
"""

#filtering
#extensions = [Extension("*", ["fe/materials/umatlibrary.pyx"],                                    
#                                 include_dirs=[numpy.get_include()] + [join(rootDirectory,umat, "include") for umat in umats], 
#                                 library_dirs= [join(rootDirectory,umat, "lib") for umat in umats] 
#                                     + [auxLibraries[lib] for lib in umatAuxLibs],
#                                 libraries= [umat for umat in umats]+[lib for lib in umatAuxLibs],    
#                                 language="c++",
#                                 extra_compile_args=compFlags)]

extensions = [Extension("*", ["fe/materials/umatlibrary.pyx"],                                    
                                 include_dirs=[numpy.get_include()] + [join(rootDirectory,'bftUserLibrary', "include")], 
                                 library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                                 runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                                 libraries= ['bftUserLibrary'],
                                 language="c++",
                                 extra_compile_args=compFlags)]

"""
Build Extensions for cython elements
"""

# filtering
cythonElements= {el : elExtraLibs for el, elExtraLibs in cythonElements.items() if exists(join(rootDirectory, el))}
for el, elementExtraLibs in cythonElements.items():
    
    # no filtering, all aux. libs must exist!
#    extraLibDirs = [auxLibraries[d] for d in elementExtraLibs] 
#    libs = [el]  + elementExtraLibs
    
#    extensions.append( Extension("*",
#                            sources=[join("fe/elements", el.lower(), "element.pyx")],
#                            language="c++",
#                            extra_compile_args=compFlags,
#                            include_dirs=[
#                                join(rootDirectory,el, "include"), numpy.get_include()],
#                            library_dirs= [join(rootDirectory,el, "lib")] + extraLibDirs,
#                            ))
    
    extensions.append( Extension("*",
                        sources=[join("fe/elements", el.lower(), "element.pyx")],
                        language="c++",
                        extra_compile_args=compFlags,
                        
                        include_dirs=[
                            join(rootDirectory,'bftUserLibrary', "include"), numpy.get_include()],
                         runtime_library_dirs= [join(rootDirectory,'bftUserLibrary', "lib") ] ,
                         libraries= ['bftUserLibrary'],
                        ))
"""
The parallel NISTP solver with OpenMP
"""
extensions += [Extension("*",
                sources = ["fe/solvers/nonlinearimplicitstaticparallel.pyx"],
                        include_dirs=[ numpy.get_include()],
                        language='c++',
                        extra_compile_args=['-fopenmp'],
                        extra_link_args=['-fopenmp'],)
                        ]  
    
"""
... And all remaining .pyx files which require no special setup
"""
extensions += [Extension("fe/*/*",
                sources = ["fe/*/*.pyx"],
                        include_dirs=[numpy.get_include()],
                        language='c++',)
                        ]
      	
"""
Compile!
"""
setup(ext_modules = cythonize(extensions,
                            compiler_directives=directives))
