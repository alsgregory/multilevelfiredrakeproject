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

Save=CopyEnsembleHierarchy(ensemble_case_example).Copy

def test_ensemble_forecast(ensemble_case_example):
    N=10
    forecast=EnsembleForecast(ensemble_case_example,ensemble_case_example.Weights,ensemble_case_example._EnsembleHierarchy__OriginalFunctionSpaces[-1][1],N)
    return forecast


Forecast=test_ensemble_forecast(ensemble_case_example)


def test_ensemble_forecast_transfer(Forecast):
    current_type=Forecast.Type
    if current_type=='Data':
        Forecast.EnsembleTransfer('Function')
        if type(Forecast.Forecast[0])!=Function:
            raise ValueError('EnsembleTransfer within EnsembleForecast doesnt work')
    if current_type=='Function':
        Forecast.EnsembleTransfer('Data')
        if type(Forecast.Forecast[0])!=np.ndarray:
            raise ValueError('EnsembleTransfer within EnsembleForecast doesnt work')

test_ensemble_forecast_transfer(Forecast)

def test_ensemble_forecast_original_ensemble(Save,Forecast):
    Save.EnsembleTransfer('Data')
    if Forecast.Ensemble[0][0][0,0]!=Save.Ensemble[0][0][0,0]:
        raise ValueError('failed')


test_ensemble_forecast_original_ensemble(Save,Forecast)


def test_ensemble_forecast_positive_definite(Forecast):
    covariance=Forecast.Cov
    if np.all(np.linalg.eigh(covariance)[0]>0)!=1:
        raise ValueError('covariance matrix is not positive definite')

test_ensemble_forecast_positive_definite(Forecast)


