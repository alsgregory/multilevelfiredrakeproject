from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level


def test_prolong_to_finest_level():
    M=MeshHierarchy(UnitSquareMesh(10,10),3)
    V=FunctionSpaceHierarchy(M,'DG',0)
    F=FunctionHierarchy(V)
    # check for prolonging to top level
    F[0].interpolate(Expression("3"))
    F[len(M)-1].interpolate(Expression("3"))
    H=FunctionHierarchy(V)
    A=ProlongInject().ProlongUpToFinestLevel(F[0],H)
    if get_level(A)[1]!=len(M)-1:
        raise AssertionError('failed')
    if norm(assemble(A-F[len(M)-1]))>0:
        raise AssertionError('hasnt prolonged. failed')

test_prolong_to_finest_level()

def test_prolong_with_non_hierarchy_function():
    M=MeshHierarchy(UnitSquareMesh(10,10),3)
    V=FunctionSpaceHierarchy(M,'DG',0)
    F=Function(V[0])
    # check for prolonging to top level
    F.interpolate(Expression("3"))
    H=FunctionHierarchy(V)
    a=0
    try:
        A=ProlongInject().ProlongUpToFinestLevel(F,H)
    except IndexError:
        a=1
    if a==0:
        raise AssertionError('failed')


test_prolong_with_non_hierarchy_function()


def test_prolong_to_any_level():
    M=MeshHierarchy(UnitSquareMesh(10,10),3)
    V=FunctionSpaceHierarchy(M,'DG',0)
    F=FunctionHierarchy(V)
    # check for prolonging to top level
    F[0].interpolate(Expression("3"))
    F[2].interpolate(Expression("3"))
    H=FunctionHierarchy(V)
    A=ProlongInject().ProlongUpToAnyLevel(2,F[0],H)
    if get_level(A)[1]!=2:
        raise AssertionError('failed')
    if norm(assemble(A-F[2]))>0:
        raise AssertionError('hasnt prolonged. failed')


test_prolong_to_any_level()


def test_inject_to_any_level():
    M=MeshHierarchy(UnitSquareMesh(10,10),3)
    V=FunctionSpaceHierarchy(M,'DG',0)
    F=FunctionHierarchy(V)
    # check for prolonging to top level
    F[-1].interpolate(Expression("3"))
    F[0].interpolate(Expression("3"))
    H=FunctionHierarchy(V)
    A=ProlongInject().InjectDownToAnyLevel(0,F[-1],H)
    if get_level(A)[1]!=0:
        raise AssertionError('failed')
    if norm(assemble(A-F[0]))>0:
        raise AssertionError('hasnt prolonged. failed')


test_inject_to_any_level()

