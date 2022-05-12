import numpy as np
import gstools
from fe.utils.misc import stringDict
from fe.utils.math import createFunction

from abc import ABC, abstractmethod


class Section(ABC):
    @abstractmethod
    def __init__(self, name, options, materialName, t, modelInfo):
        pass

    @abstractmethod
    def assignSectionPropertiesToModel(self, modelInfo):
        pass
