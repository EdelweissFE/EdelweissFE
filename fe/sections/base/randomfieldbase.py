import numpy as np
import gstools
from fe.utils.misc import stringDict
from fe.utils.math import createFunction
from fe.sections.base.sectionbase import Section as SectionBase

from abc import ABC, abstractmethod


class Section(SectionBase):
    def __init__(self, name, options, materialName, t, modelInfo):

        self.options = stringDict([e for l in options for e in l])
        options = self.options

        self.materialName = materialName

        self.elSetNames = [s.strip() for s in options["elSets"].split(",")]

        dimension = modelInfo["domainSize"]
        variance = float(options["variance"])
        lengthScale = float(options["lengthScale"])
        seed = int(options["seed"])

        model = gstools.Gaussian(
            dim=dimension,
            var=variance,
            len_scale=lengthScale,
        )
        self.srf = gstools.SRF(model, seed=seed)

        self.randomFunction = createFunction(options["f(x,ref,rand)"], "x", "ref", "rand", modelInfo=modelInfo)
