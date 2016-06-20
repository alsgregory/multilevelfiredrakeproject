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


def test_ensemble_transfer_consistency():

    current_type = ensemble_case_example.Type

    if current_type == 'Data':
        Save = CopyEnsembleHierarchy(ensemble_case_example).Copy

    else:

        assert current_type == 'Function'

        ensemble_case_example.EnsembleTransfer('Data')
        Save = CopyEnsembleHierarchy(ensemble_case_example).Copy

    # Convert twice
    ensemble_case_example.EnsembleTransfer('Function')

    assert ensemble_case_example.Type == 'Function'

    ensemble_case_example.EnsembleTransfer('Data')

    assert ensemble_case_example.Type == 'Data'

    assert ensemble_case_example.Ensemble[0][0][0, 0] == Save.Ensemble[0][0][0, 0]


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
