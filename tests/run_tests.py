''' Run Tests '''


from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
import pytest


# Now run tests

import test_import

import test_prolong_inject

import test_ensemble_forecast

import test_ensemble_transfer

import test_ensemble_hierarchy

import test_sample_statistics

import test_state

import test_convergence







