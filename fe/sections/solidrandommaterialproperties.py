import numpy as np
from fe.sections.base.randomfieldbase import Section as RandomBase


class Section(RandomBase):
    def __init__(self, name, options, materialName, t, model):
        super().__init__(name, options, materialName, t, model)

        self.indexRandom = int(self.options["indexRandom"])

    def assignSectionPropertiesToModel(self, model):
        elSets = [model.elementSets[setName] for setName in self.elSetNames]
        material = model.materials[self.materialName]

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

        return model
