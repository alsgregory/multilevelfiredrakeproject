from __future__ import division  # Get proper divison
import numpy as np
import random


from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from test_ensemble_case import *


def test_ensemble_hierarchy_generation():

    ensemble_hierarchy = EnsembleHierarchy()

    assert ensemble_hierarchy.Ensemble == []


def test_ensemble_appending():

    n = 100
    L = 1

    mesh = UnitSquareMesh(5, 5)
    Mesh_Hierarchy = GenerateMeshHierarchy(mesh, L)

    V = FunctionSpaceHierarchy(Mesh_Hierarchy, 'DG', 0)

    assert len(Mesh_Hierarchy) == L + 2

    ensemble_hierarchy = EnsembleHierarchy()

    l = 0
    for i in range(n):

        u = FunctionHierarchy(V)
        v = FunctionHierarchy(V)

        A = tuple([u[-1], v[-1]])  # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A, l)

    l = 1
    for i in range(1):

        u = FunctionHierarchy(V)
        v = FunctionHierarchy(V)

        A = tuple([u[-1], v[-1]])  # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A, l)

    assert len(ensemble_hierarchy.Ensemble) == 2

    assert len(ensemble_hierarchy.Ensemble[0]) == n

    assert len(ensemble_hierarchy.Ensemble[1]) == 1


def test_ensemble_same_function_space():

    n = 1
    L = 1

    mesh = UnitSquareMesh(5, 5)
    Mesh_Hierarchy = GenerateMeshHierarchy(mesh, L)

    V = FunctionSpaceHierarchy(Mesh_Hierarchy, 'DG', 0)

    ensemble_hierarchy = EnsembleHierarchy()

    l = 0
    a = 0
    for i in range(n):

        u = FunctionHierarchy(V)
        v = FunctionHierarchy(V)

        A = tuple([u[-2], v[-1]])  # tuple of same function space

        try:
            ensemble_hierarchy.AppendToEnsemble(A, l)

        except:
            a = 1

    assert a != 0


def test_ensemble_type():

    n = 1
    L = 1

    mesh = UnitSquareMesh(5, 5)
    Mesh_Hierarchy = GenerateMeshHierarchy(mesh, L)

    V = FunctionSpaceHierarchy(Mesh_Hierarchy, 'DG', 0)

    ensemble_hierarchy = EnsembleHierarchy()

    l = 0
    for i in range(n):

        u = FunctionHierarchy(V)
        v = FunctionHierarchy(V)

        A = tuple([u[-1], v[-1]])
        ensemble_hierarchy.AppendToEnsemble(A, l)

    ensemble_hierarchy.EnsembleTransfer('Data')

    assert np.shape(ensemble_hierarchy.Ensemble) == (
        1, 2, len(u[-1].dat.data), n)

    assert ensemble_hierarchy.Type == 'Data'


def test_ensemble_nxl():

    n = 1
    L = 1

    mesh = UnitSquareMesh(5, 5)
    Mesh_Hierarchy = GenerateMeshHierarchy(mesh, L)

    V = FunctionSpaceHierarchy(Mesh_Hierarchy, 'DG', 0)

    assert len(Mesh_Hierarchy) == L + 2

    ensemble_hierarchy = EnsembleHierarchy()

    l = 0
    for i in range(n):

        u = FunctionHierarchy(V)
        v = FunctionHierarchy(V)

        A = tuple([u[-1], v[-1]])  # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A, l)

    l = 1
    for i in range(n):

        u = FunctionHierarchy(V)
        v = FunctionHierarchy(V)

        A = tuple([u[-1], v[-1]])  # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A, l)

    # test nxl
    assert hasattr(ensemble_hierarchy, 'nxl') == 1

    assert len(ensemble_hierarchy.nxl) == 2

    assert ensemble_hierarchy.nxl[0] == (5 * 5)

    assert ensemble_hierarchy.nxl[1] == (10 * 10)


if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
