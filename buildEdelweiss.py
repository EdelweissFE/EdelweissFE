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
              'wraparound':False}

gccFlags = ["-std=c++11", "-w", "-g3"]
msvcFlags = []

compFlags = gccFlags

# rootDirectory where the umats are located
rootDirectory = expanduser("~/Dropbox/PHD/GITrepositoryLinux")

# 0.)
# list of directories, where the extra libs to link are located (e.g. the directory
# of bftMechancics)
# + the list of extra libraries to link (e.g. bftMechanics)

extraLibs = [#(join(rootDirectory, 'SolidificationCreep/lib'), 'SolidificationCreep'),
            #(join(rootDirectory, 'ModLeonCore/lib'),'ModLeonCore'),
            (join(rootDirectory, 'bftMechanics/lib'),'bftMechanics'),
            # (join(rootDirectory, 'bftMaterialLibrary/lib'),'bftMaterialLibrary'),
            ]

# we check if every directory exists ... or else we don't include it
extraLibs = [libDef for libDef in extraLibs if exists(libDef[0])]  
umats = [
#         "Meschke",
         "ModLeon",
        "ModLeonNonLocal",
        # "ModLeonPlaneStress",
        # "ModLeonPlaneStressWrapped",
        # "ModLeonAdaptive",
        # "ModLeonAnalytical",
        # "ModLeonSemiExplicit",
        # "ModLeonSemiExplicitAdaptive",
        # "ShotLeon",
        # "SchaedlichSchweiger",
        # "Schuetz",
        # "linearElasticSolidificationCreep",
        # "Unteregger",
        ]
        
# we check if every directory exists ... or else we don't include it
umats = [umat for umat in umats if exists(join(rootDirectory, umat))]        
extensions = []            
#                        
extensions += [Extension("*", ["fe/materials/umatlibrary.pyx"],                                    
                                 include_dirs=[numpy.get_include()] + [join(rootDirectory,umat, "include") for umat in umats], 
                                 library_dirs=                           # linkDirlist = extralibs list + the UMAT.lib directory
                                 [join(rootDirectory,umat, "lib") for umat in umats] + [d for d, lib in extraLibs],
                                 libraries= [umat for umat in umats]+[lib for d, lib in extraLibs],    
                                 language="c++",
                                 extra_compile_args=compFlags)]
# 1.b.)
# create a UMAT-extension definition for each UMAT in list.

cythonElements = [  'uelCPE4',
                    # 'uelCPE4R',
                    'uelCPS4', 
                    'uelCPS4nonLocal',
                    
                    # 'uelCPS8R', 
                    'uelCPS8RNonLocal', 
                    ]

cythonElements= [el for el in cythonElements if exists(join(rootDirectory, el))]        

extensions += [Extension("*",
                        sources=[join("fe/elements", el.lower(), "element.pyx")],
                        language="c++",
                        extra_compile_args=compFlags,
                        include_dirs=[
                            join(rootDirectory,el, "include"),                  # list of all include directories, including UMAT headers
                                      numpy.get_include()],
                        library_dirs= [d for d, lib in extraLibs]+                          # linkDirlist = extralibs list + the UMAT.lib directory
                                      [join(rootDirectory,el, "lib")],
                        libraries= [el]+[lib for d, lib in extraLibs]+[lib for d, lib in extraLibs] 
                        )                    # linkList = extralibs + UMAT.lib  
                        for el in cythonElements]                   
      	
# 3) compile!
setup(ext_modules = cythonize(extensions,
                            compiler_directives=directives))
