import numpy as np
import gstools
from fe.utils.misc import stringDict


class Section:
    def __init__(self, name, options, materialName, t, modelInfo):
        self.elSetNames = [e for l in options for e in l]
        self.materialName = materialName

    def assignSectionPropertiesToModel(self, modelInfo):
        elSets = [modelInfo["elementSets"][setName] for setName in self.elSetNames]
        material = modelInfo["materials"][self.materialName]

        for elSet in elSets:
            for el in elSet:
                el.initializeElement()
                el.setMaterial(material["name"], material["properties"])

        return modelInfo
