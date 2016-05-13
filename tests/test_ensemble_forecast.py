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
    forecast=EnsembleForecast(ensemble_case_example,ensemble_case_example._EnsembleHierarchy__OriginalFunctionSpaces[-1][1])
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

def test_ensemble_forecast_delta(Forecast):
    if Forecast.Delta>1 or Forecast.Delta<0:
        raise ValueError('Delta used to shrink covariance matrix is not within [0,1]')

test_ensemble_forecast_delta(Forecast)


def test_ensemble_forecast_Flatten_Array(Forecast):
    A=np.ones((10,10,10))
    B=Forecast._EnsembleForecast__FlattenMultiD(A)
    if np.shape(B)!=(100,10):
        raise ValueError('Flatten Array functionality fails')
    return B


B=test_ensemble_forecast_Flatten_Array(Forecast)


def test_ensemble_forecast_Return_Array(Forecast,B):
    A=Forecast._EnsembleForecast__ReturnMultiD(B,(10,10))
    if np.shape(A)!=(10,10,10):
        raise ValueError('Return Multi Dimension functionality fails')


test_ensemble_forecast_Return_Array(Forecast,B)

