from fe.utils.misc import convertLinesToStringDictionary
from fe.config.solvers import getSolverByName
from fe.journal.journal import Journal


def createSolversFromInputFile(inputfile: dict, jobInfo: dict, journal: Journal):
    solvers = {}
    for solverDefinition in inputfile["*solver"]:
        solverType = solverDefinition["solver"]
        solverName = solverDefinition["name"]
        solverData = solverDefinition["data"]

        Solver = getSolverByName(solverType)

        solverData = convertLinesToStringDictionary(solverData)

        solvers[solverName] = Solver(jobInfo, journal, **solverData)

    return solvers
