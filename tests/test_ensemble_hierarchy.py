from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from test_ensemble_case import *

ensemble_case_example = make_ensemble()


def test_ensemble_hierarchy_generation():
    ensemble_hierarchy=EnsembleHierarchy()
    if ensemble_hierarchy.Ensemble!=[]:
        raise ValueError('failed')

test_ensemble_hierarchy_generation()



def test_ensemble_appending():
    n=100
    l=0
    ensemble_hierarchy=EnsembleHierarchy()
    for i in range(n):
        A=[Function(FunctionSpace(UnitSquareMesh(5,5),'DG',0)).assign(0),Function(FunctionSpace(UnitSquareMesh(10,10),'DG',0)).assign(0)]
        ensemble_hierarchy.AppendToEnsemble(A,l)
    l=1
    for i in range(1):
        A=[Function(FunctionSpace(UnitSquareMesh(5,5),'DG',0)).assign(0),Function(FunctionSpace(UnitSquareMesh(10,10),'DG',0)).assign(0)]
        ensemble_hierarchy.AppendToEnsemble(A,l)
    if len(ensemble_hierarchy.Ensemble)!=2:
        raise ValueError('failed')
    if len(ensemble_hierarchy.Ensemble[0])!=100:
        raise ValueError('failed')
    if len(ensemble_hierarchy.Ensemble[1])!=1:
        raise ValueError('failed')


test_ensemble_appending()


def test_ensemble_type():
    n=100
    l=0
    ensemble_hierarchy=EnsembleHierarchy()
    for i in range(n):
        A=[Function(FunctionSpace(UnitSquareMesh(10,10),'DG',0)).assign(0),Function(FunctionSpace(UnitSquareMesh(10,10),'DG',0)).assign(0)]
        ensemble_hierarchy.AppendToEnsemble(A,l)
    E=EnsembleHierarchyTransfer(ensemble_hierarchy)
    E.EnsembleTransfer()
    if np.shape(ensemble_hierarchy.Ensemble)!=(1,2,len(Function(FunctionSpace(UnitSquareMesh(10,10),'DG',0)).dat.data),n):
        raise ValueError('failed')
    if ensemble_hierarchy.Type!='Data':
        raise ValueError('failed')


test_ensemble_type()





