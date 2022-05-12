import numpy as np
import gstools
from fe.utils.misc import stringDict
from fe.utils.math import createFunction


class Section:
    def __init__(self, name, options, materialName, t, modelInfo):

        self.options = stringDict([e for l in options for e in l])
        options = self.options

        self.materialName = materialName

        self.elSetNames = [options["elSet"]]
        self.referenceThickness = np.array([float(options["thickness"])], dtype=float)

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

    def assignSectionPropertiesToModel(self, modelInfo):

        elSets = [modelInfo["elementSets"][setName] for setName in self.elSetNames]
        material = modelInfo["materials"][self.materialName]

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

        return modelInfo
