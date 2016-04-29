from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level

from test_problem_functions import *

from test_ensemble_case import *

ensemble_case_example = make_ensemble()

def test_ensemble_transfer_consistency(ensemble_case_example):
    current_type=ensemble_case_example.Type
    if current_type=='Data':
        Save=CopyEnsembleHierarchy(ensemble_case_example).Copy
    else:
        if current_type!='Function':
            raise AttributeError('failed. type doesnt exist')
        else:
            ensemble_case_example.EnsembleTransfer('Data')
            Save=CopyEnsembleHierarchy(ensemble_case_example).Copy
    # Convert twice
    ensemble_case_example.EnsembleTransfer('Function')
    if ensemble_case_example.Type!='Function':
        raise ValueError('failed. type hasnt changed')
    else:
        ensemble_case_example.EnsembleTransfer('Data')
        if ensemble_case_example.Type!='Data':
            raise ValueError('failed. type hasnt changed')
        else:
            if ensemble_case_example.Ensemble[0][0][0,0]!=Save.Ensemble[0][0][0,0]:
                raise ValueError('failed')


test_ensemble_transfer_consistency(ensemble_case_example)



