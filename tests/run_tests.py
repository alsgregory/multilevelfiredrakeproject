''' Run Tests '''


from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *


# Now run tests

import test_import
import test_discretization # EDIT!
import test_ensemble_transfer
import test_sample_statistics


#import test_ensemble_hierarchy # EDIT!





