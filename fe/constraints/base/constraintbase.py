from abc import ABC, abstractmethod
from fe.variables.scalarvariable import ScalarVariable
from numpy import ndarray


class ConstraintBase(ABC):
    """The constraint base class.

    Parameters
    ----------
    name
        The name of the constraint.
    options
        A dictionary containing the options for the constraint.
    modelInfo
        A dictionary containing the model tree.
    """

    @abstractmethod
    def __init__(self, name: str, options: dict, modelInfo: dict):

        self.name = name
        self.nodes = []
        self.fieldsOnNodes = [
            [],
        ]
        self.scalarVariables = []
        self.nDof = 0

    @abstractmethod
    def getNumberOfAdditionalNeededScalarVariables(
        self,
    ) -> int:
        """This method is called to determine the additional number of scalar variables
        this Constraint needs.

        Returns
        -------
        int
            Number of ScalarVariables required.

        """

        pass

    def assignAdditionalScalarVariables(self, scalarVariables: list[ScalarVariable]):
        """Assign a list of scalar variables associated with this constraint.

        Parameters
        ----------
        list
            The list of ScalarVariables associated with this constraint.

        """

        self.scalarVariables = scalarVariables

    @abstractmethod
    def applyConstraint(self, Un1: ndarray, dU: ndarray, PExt: ndarray, V: ndarray, increment: tuple):
        """Apply the constraint.  Add the contributions to the external load vector and the system matrix.

        Parameters
        ----------
        Un1
            The current total solution vector.
        dU
            The current increment since the last time the constraint was applied.
        PExt
            The external load vector.
        V
            The system matrix in vector form.
        increment
            The time increment information.
        """

        pass
