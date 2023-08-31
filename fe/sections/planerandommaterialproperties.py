import numpy as np
import gstools
from fe.sections.base.randomfieldbase import Section as RandomBase


class Section(RandomBase):
    def __init__(self, name, options, materialName, t, model):
        super().__init__(name, options, materialName, t, model)

        self.thickness = np.array(
            [
                float((self.options["thickness"])),
            ]
        )
        self.indexRandom = int(self.options["indexRandom"])

    def assignSectionPropertiesToModel(self, model):
        elSets = [model.elementSets[setName] for setName in self.elSetNames]
        material = model.materials[self.materialName]

        elementThickness = self.thickness

        indexRandom = self.indexRandom

        for elSet in elSets:
            for el in elSet:
                elProperties = elementThickness
                el.setProperties(elProperties)
                el.initializeElement()

                x = el.getCoordinatesAtCenter()
                randomFieldValue = self.srf(x)

                materialProperties = np.copy(material["properties"])
                materialProperties[indexRandom] = self.randomFunction(
                    x, materialProperties[indexRandom], randomFieldValue
                )

                el.setMaterial(material["name"], materialProperties)

        return model
