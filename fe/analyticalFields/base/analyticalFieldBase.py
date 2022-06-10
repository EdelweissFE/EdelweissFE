# import numpy as np
# import gstools
# from fe.utils.misc import stringDict
# from fe.utils.math import createFunction
#
from abc import ABC, abstractmethod


class AnalyticalField(ABC):
    @abstractmethod
    def __init__(self, name, data, modelInfo):
        pass
