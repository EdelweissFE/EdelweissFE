from abc import ABC, abstractmethod
from fe.elements.scalarvariable import ScalarVariable


class ConstraintBase(ABC):
    @abstractmethod
    def __init__(self, name, options, modelInfo):
        self.name = name
        self.nodes = []
        self.fieldsOnNodes = [
            [],
        ]
        self.scalarVariables = []
        self.nDof = 0

    # @abstractmethod
    # def update(self, options : dict[str], modelInfo : dict):
    #     """Update the options of a constraint.

    #     Parameters
    #     ----------
    #     options
    #         options in string format.
    #     modelInfo
    #         The modelInfo tree.

    #     """

    #     pass

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

    @abstractmethod
    def assignAdditionalScalarVariables(self, scalarVariables: list[ScalarVariable]):
        """Assign a list of scalar variables associated with this constraint.

        Parameters
        ----------
        list
            The list of ScalarVariables associated with this constraint.

        """
        self.scalarVariables = scalarVariables
