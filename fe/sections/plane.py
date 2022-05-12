import numpy as np
import gstools
from fe.utils.misc import stringDict


class Section:
    def __init__(self, name, options, materialName, thickness, modelInfo):

        self.materialName = materialName
        self.elSetNames = [e for l in options for e in l]
        self.thickness = thickness

    def assignSectionPropertiesToModel(self, modelInfo):

        elSets = [modelInfo["elementSets"][setName] for setName in self.elSetNames]
        material = modelInfo["materials"][self.materialName]

        for elSet in elSets:
            for el in elSet:
                elementThickness = self.thickness
                elProperties = np.array(
                    [
                        elementThickness,
                    ],
                    dtype=np.float,
                )
                el.setProperties(elProperties)
                el.initializeElement()
                el.setMaterial(material["name"], material["properties"])

        return modelInfo
