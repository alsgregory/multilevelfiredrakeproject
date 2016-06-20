from __future__ import division  # Get proper divison
import numpy as np
import random
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *


def test_state_levels_1():

    M = UnitSquareMesh(10, 10)
    MH = MeshHierarchy(M, 1)

    V = FunctionSpaceHierarchy(MH, 'DG', 0)

    F = FunctionHierarchy(V)

    b = 1
    a = 0

    try:
        Sol = state(F[0], F[1])

    except Warning:
        a = 1
        b = 0

    assert a < b


def test_state_levels_2():

    M = UnitSquareMesh(10, 10)
    V = FunctionSpace(M, 'DG', 0)

    F = Function(V)
    G = Function(V)

    b = 1
    a = 0

    try:
        Sol = state(F, G)

    except Warning:
        a = 1
        b = 0

    assert a > b


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
