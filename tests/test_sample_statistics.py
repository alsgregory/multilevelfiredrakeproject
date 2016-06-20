from __future__ import division # Get proper divison
import numpy as np
import random

from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level

from test_problem_functions import *

from test_ensemble_case import *

ensemble_case_example = make_ensemble()

def test_sample_statistics_length():
    SampleStatistics(ensemble_case_example)
    # check for length
    L=ensemble_case_example.L
    assert len(ensemble_case_example.Mean)==L
    assert len(ensemble_case_example.Variance)==L

ensemble_case_example = make_ensemble()

def test_sample_statistics_type():
    ensemble_case_example.EnsembleTransfer('Data')
    SampleStatistics(ensemble_case_example)
    a=0
    if ensemble_case_example.Type=='Data':
        assert type(ensemble_case_example.MultilevelExpectation)==Function
        assert type(ensemble_case_example.Mean[0])==Function
        assert type(ensemble_case_example.Variance[0])==Function
    else:
        a=1
        assert a==0
    ensemble_case_example.EnsembleTransfer('Function')
    SampleStatistics(ensemble_case_example)
    if ensemble_case_example.Type=='Function':
        assert type(ensemble_case_example.MultilevelExpectation)==Function
        assert type(ensemble_case_example.Mean[0])==Function
        assert type(ensemble_case_example.Variance[0])==Function
    else:
        a=1
        assert a==0

if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
