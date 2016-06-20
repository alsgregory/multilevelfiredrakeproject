from __future__ import division # Get proper divison
import numpy as np
import random

from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
import sys
from firedrake.mg.utils import get_level


from test_problem_functions import *


'''  Test Imports '''


def module_check():
    # test the imports
    assert hasattr(EnsembleHierarchy,'EnsembleTransfer')==True
    assert hasattr(EnsembleHierarchy,'AppendToEnsemble')==True
    assert hasattr(EnsembleForecast,'EnsembleTransfer')==True

if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))

