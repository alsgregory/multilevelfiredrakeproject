from __future__ import division  # Get proper divison
import numpy as np
import random

from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level

from test_problem_functions import *

from test_ensemble_case import *


ensemble_case_example = make_ensemble()

Save = CopyEnsembleHierarchy(ensemble_case_example).Copy


def test_ensemble_forecast():

    forecast = EnsembleForecast(
        ensemble_case_example, ensemble_case_example._EnsembleHierarchy__OriginalFunctionSpaces[-1][1])

    return forecast


Forecast = test_ensemble_forecast()


def test_ensemble_forecast_transfer():

    current_type = Forecast.Type

    if current_type == 'Data':
        Forecast.EnsembleTransfer('Function')

        assert isinstance(Forecast.Forecast[0], Function)

    if current_type == 'Function':
        Forecast.EnsembleTransfer('Data')

        assert isinstance(Forecast.Forecast[0], np.ndarray)


def test_ensemble_forecast_original_ensemble():

    Save.EnsembleTransfer('Data')

    assert Forecast.Ensemble[0][0][0, 0] == Save.Ensemble[0][0][0, 0]


def test_ensemble_forecast_positive_definite():

    covariance = Forecast.Cov

    assert np.all(np.linalg.eigh(covariance)[0] > 0) == 1


def test_ensemble_forecast_diag_positivity():

    covariance = Forecast.Cov

    assert np.all(np.diag(covariance) > 0) == 1


def test_ensemble_forecast_delta():

    assert Forecast.Delta <= 1

    assert Forecast.Delta >= 0


def test_ensemble_forecast_Flatten_Array():

    A = np.ones((10, 10, 10))
    B = Forecast._EnsembleForecast__FlattenMultiD(A)

    assert np.shape(B) == (100, 10)

    return B


B = test_ensemble_forecast_Flatten_Array()


def test_ensemble_forecast_Return_Array():

    A = Forecast._EnsembleForecast__ReturnMultiD(B, (10, 10))

    assert np.shape(A) == (10, 10, 10)

if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
