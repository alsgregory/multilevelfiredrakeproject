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

def test_sample_statistics_length(ensemble_case_example):
    SampleStatistics(ensemble_case_example)
    # check for length
    L=ensemble_case_example.L
    if len(ensemble_case_example.Mean)!=L:
        raise ValueError('failed')
    if len(ensemble_case_example.Variance)!=L:
        raise ValueError('failed')


test_sample_statistics_length(ensemble_case_example)



def test_sample_statistics_type(ensemble_case_example):
    if ensemble_case_example.Type=='Data':
        if type(ensemble_case_example.MultilevelExpectation)!=np.ndarray:
            raise ValueError('failed')
        if type(ensemble_case_example.Mean[0])!=np.ndarray:
            raise ValueError('failed')
        if type(ensemble_case_example.Variance[0])!=np.ndarray:
            raise ValueError('failed')
    if ensemble_case_example.Type=='Function':
        if type(ensemble_case_example.MultilevelExpectation)!=Function:
            raise ValueError('failed')
        if type(ensemble_case_example.Mean[0])!=Function:
            raise ValueError('failed')
        if type(ensemble_case_example.Variance[0])!=Function:
            raise ValueError('failed')

test_sample_statistics_type(ensemble_case_example)

