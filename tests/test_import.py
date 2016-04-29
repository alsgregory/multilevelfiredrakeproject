from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
import sys
from firedrake.mg.utils import get_level


from test_problem_functions import *


'''  Test Imports '''


def module_check():
    if hasattr(EnsembleHierarchy,'EnsembleTransfer')==False:
        raise AttributeError('failed')
    if hasattr(EnsembleHierarchy,'AppendToEnsemble')==False:
        raise AttributeError('failed')
    #
    if hasattr(EnsembleForecast,'EnsembleTransfer')==False:
        raise AttributeError('failed')


module_check()
