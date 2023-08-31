from abc import ABC, abstractmethod
from fe.variables.scalarvariable import ScalarVariable
from fe.models.femodel import FEModel
from numpy import ndarray


class ConstraintBase(ABC):
    @abstractmethod
    def __init__(self, name: str, options: dict, model: FEModel):
        """The constraint base class.

        Constraints can act on nodal variables, and scalar variables.
        If scalar variables are required, the can be created on demand by
        defining
        :func:`~ConstraintBase.getNumberOfAdditionalNeededScalarVariables` and
        :func:`~ConstraintBase.assignAdditionalScalarVariables`, which are called at the beginning of an analysis.

        If scalar variables are used, EdelweissFE expects the layout of the external load vector PExt (and the stiffness)
        to be of the form

        .. code-block:: console

            [ node 1 - dofs field 1,
              node 1 - dofs field 2,
              node 1 - ... ,
              node 1 - dofs field n,
              node 2 - dofs field 1,
              ... ,
              node N - dofs field n,
              scalar variable 1,
              scalar variable 2,
              ... ,
              scalar variable J ].

        Parameters
        ----------
        name
            The name of the constraint.
        options
            A dictionary containing the options for the constraint.
        model
            A dictionary containing the model tree.
        """

        self.scalarVariables = []

    @property
    @abstractmethod
    def nodes(self) -> list:
        """The nodes this constraint is acting on."""

        pass

    @property
    @abstractmethod
    def fieldsOnNodes(self) -> list:
        """The fields on the nodes this constraint is acting on."""

        pass

    @property
    @abstractmethod
    def nDof(self) -> int:
        """The total number of degrees of freedom this constraint is associated with."""

        pass

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

        return 0

    def assignAdditionalScalarVariables(self, scalarVariables: list[ScalarVariable]):
        """Assign a list of scalar variables associated with this constraint.

        Parameters
        ----------
        list
            The list of ScalarVariables associated with this constraint.

        """

        self.scalarVariables = scalarVariables

    @abstractmethod
    def applyConstraint(self, U_np: ndarray, dU: ndarray, PExt: ndarray, V: ndarray, increment: tuple):
        """Apply the constraint.  Add the contributions to the external load vector and the system matrix.

        Parameters
        ----------
        U_np
            The current total solution vector.
        dU
            The current increment since the last time the constraint was applied.
        PExt
            The external load vector.
        K
            The system (stiffness) matrix.
        increment
            The time increment information.
        """

        pass
