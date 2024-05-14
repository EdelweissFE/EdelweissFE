import gstools
import numpy as np

from edelweissfe.utils.math import createFunction
from edelweissfe.utils.misc import convertLinesToStringDictionary


class Section:
    def __init__(self, name, options, materialName, t, model):
        self.options = convertLinesToStringDictionary(options)
        options = self.options

        self.materialName = materialName

        self.elSetNames = [options["elSet"]]
        self.referenceThickness = np.array([float(options["thickness"])], dtype=float)

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

    def assignSectionPropertiesToModel(self, model):
        elSets = [model.elementSets[setName] for setName in self.elSetNames]
        material = model.materials[self.materialName]

        for elSet in elSets:
            for el in elSet:
                # elementThickness = self.referenceThickness * (1 + self.srf(el.getCoordinatesAtCenter()))

                x = el.getCoordinatesAtCenter()
                randomFieldValue = self.srf(x)
                elementThickness = self.randomFunction(x, self.referenceThickness, randomFieldValue)

                elProperties = elementThickness
                el.setProperties(elProperties)
                el.initializeElement()
                el.setMaterial(material["name"], material["properties"])

        return model
