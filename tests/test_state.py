from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *

def test_state_levels_1():
    M=UnitSquareMesh(10,10)
    MH=MeshHierarchy(M,1)
    V=FunctionSpaceHierarchy(MH,'DG',0)
    F=FunctionHierarchy(V)
    try:
        Sol=state(F[0],F[1])
    except Warning:
        raise ValueError('failed')

test_state_levels_1()

def test_state_levels_2():
    M=UnitSquareMesh(10,10)
    V=FunctionSpace(M,'DG',0)
    F=Function(V)
    G=Function(V)
    b=1; a=0
    try:
        Sol=state(F,G)
    except Warning:
        a=1; b=0
    if a<b:
        raise ValueError('failed')

test_state_levels_2()


