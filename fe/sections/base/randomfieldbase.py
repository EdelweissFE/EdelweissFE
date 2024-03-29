import numpy as np
import gstools
from fe.utils.misc import convertLinesToStringDictionary
from fe.utils.math import createFunction
from fe.sections.base.sectionbase import Section as SectionBase

import warnings


class Section(SectionBase):
    def __init__(self, name, options, materialName, t, model):
        warnings.warn(
            "'solidrandommaterialproperties' section type is deprecated and will be removed in the future; use 'solid' or 'plane' section section with 'materialParametersFromField' option instead"
        )

        self.options = convertLinesToStringDictionary(options)
        options = self.options

        self.materialName = materialName

        self.elSetNames = [s.strip() for s in options["elSets"].split(",")]

        dimension = model.domainSize
        variance = float(options["variance"])
        lengthScale = float(options["lengthScale"])
        seed = int(options["seed"])

        model = gstools.Gaussian(
            dim=dimension,
            var=variance,
            len_scale=lengthScale,
        )
        self.srf = gstools.SRF(model, seed=seed)

        self.randomFunction = createFunction(options["f(x,ref,rand)"], "x", "ref", "rand", model=model)
