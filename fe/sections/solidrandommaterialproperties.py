import numpy as np
import gstools
from fe.utils.misc import stringDict
from fe.sections.base.randomfieldbase import Section as RandomBase


class Section(RandomBase):
    def __init__(self, name, options, materialName, t, modelInfo):
        super().__init__(name, options, materialName, t, modelInfo)

        self.indexRandom = int(self.options["indexRandom"])

    def assignSectionPropertiesToModel(self, modelInfo):

        elSets = [modelInfo["elementSets"][setName] for setName in self.elSetNames]
        material = modelInfo["materials"][self.materialName]

        indexRandom = self.indexRandom

        for elSet in elSets:
            for el in elSet:
                el.initializeElement()

                x = el.getCoordinatesAtCenter()
                randomFieldValue = self.srf(x)

                materialProperties = np.copy(material["properties"])
                materialProperties[indexRandom] = self.randomFunction(
                    x, materialProperties[indexRandom], randomFieldValue
                )

                el.setMaterial(material["name"], materialProperties)

        return modelInfo
