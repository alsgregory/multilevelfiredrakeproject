''' Run Tests '''


from __future__ import division  # Get proper divison
import numpy as np
import random

from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *


import os

os.system('python test_import.py')

os.system('python test_prolong_inject.py')

os.system('python test_ensemble_forecast.py')

os.system('python test_ensemble_transfer.py')

os.system('python test_ensemble_hierarchy.py')

os.system('python test_sample_statistics.py')

os.system('python test_state.py')

os.system('python test_convergence.py')
