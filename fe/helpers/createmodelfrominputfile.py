from fe.models.femodel import FEModel
from fe.journal.journal import Journal
from fe.config.generators import getGeneratorFunction
from fe.utils.abqmodelconstructor import AbqModelConstructor
from fe.variables.scalarvariable import ScalarVariable


def fillFEModelFromInputFile(model: FEModel, inputfile: dict, journal: Journal):
    """Convenience helper function
    to fill an existing (possibly empty) FEModel using the input file and generators.

    Parameters
    ----------
    FEModel
        The model tree to be filled.
    input
        The processed inputfile in dictionary form.
    journal
        The Journal for logging purposes.

    Returns
    -------
    FEModel
        The updated, filled model tree.
    """

    # call individual optional model generators
    for generatorDefinition in inputfile["*modelGenerator"]:
        if generatorDefinition.get("executeAfterManualGeneration", False):
            continue
        gen = generatorDefinition["generator"]
        model = getGeneratorFunction(gen)(generatorDefinition, model, journal)

    # the standard 'Abaqus like' model generator is invoked unconditionally, and it has direct access to the inputfile
    abqModelConstructor = AbqModelConstructor(journal)
    model = abqModelConstructor.createGeometryFromInputFile(model, inputfile)
    model = abqModelConstructor.createMaterialsFromInputFile(model, inputfile)
    model = abqModelConstructor.createConstraintsFromInputFile(model, inputfile)
    model = abqModelConstructor.createAnalyticalFieldsFromInputFile(model, inputfile)
    model = abqModelConstructor.createSectionsFromInputFile(model, inputfile)

    # call individual optional model generators,
    for generatorDefinition in inputfile["*modelGenerator"]:
        if not generatorDefinition.get("executeAfterManualGeneration", False):
            continue
        gen = generatorDefinition["generator"]
        model = getGeneratorFunction(gen)(generatorDefinition, model, journal)

    return model
