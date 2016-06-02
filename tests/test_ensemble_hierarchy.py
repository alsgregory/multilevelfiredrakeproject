from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from test_ensemble_case import *


def test_ensemble_hierarchy_generation():
    ensemble_hierarchy=EnsembleHierarchy()
    if ensemble_hierarchy.Ensemble!=[]:
        raise ValueError('failed')

test_ensemble_hierarchy_generation()


def test_ensemble_appending():
    n=100
    L=1
    mesh=UnitSquareMesh(5,5)
    Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
    V=FunctionSpaceHierarchy(Mesh_Hierarchy,'DG',0)
    #
    if len(Mesh_Hierarchy)!=L+2:
        raise AssertionError('Mesh Hierarchy isnt of length L+2')
    #
    ensemble_hierarchy=EnsembleHierarchy()
    l=0
    for i in range(n):
        u=FunctionHierarchy(V)
        v=FunctionHierarchy(V)
        A=tuple([u[-1],v[-1]]) # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A,l)
    l=1
    for i in range(1):
        u=FunctionHierarchy(V)
        v=FunctionHierarchy(V)
        A=tuple([u[-1],v[-1]]) # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A,l)
    if len(ensemble_hierarchy.Ensemble)!=2:
        raise ValueError('failed')
    if len(ensemble_hierarchy.Ensemble[0])!=n:
        raise ValueError('failed')
    if len(ensemble_hierarchy.Ensemble[1])!=1:
        raise ValueError('failed')


test_ensemble_appending()



def test_ensemble_same_function_space():
    n=1
    L=1
    mesh=UnitSquareMesh(5,5)
    Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
    V=FunctionSpaceHierarchy(Mesh_Hierarchy,'DG',0)
    ensemble_hierarchy=EnsembleHierarchy()
    l=0
    a=0
    for i in range(n):
        u=FunctionHierarchy(V)
        v=FunctionHierarchy(V)
        A=tuple([u[-2],v[-1]]) # tuple of same function space
        try:
            ensemble_hierarchy.AppendToEnsemble(A,l)
        except:
            a=1
    if a==0:
        raise AssertionError('failed. same function space error hasnt arrised')
        
        


test_ensemble_same_function_space()




def test_ensemble_type():
    n=1
    L=1
    mesh=UnitSquareMesh(5,5)
    Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
    V=FunctionSpaceHierarchy(Mesh_Hierarchy,'DG',0)
    ensemble_hierarchy=EnsembleHierarchy()
    l=0
    for i in range(n):
        u=FunctionHierarchy(V)
        v=FunctionHierarchy(V)
        A=tuple([u[-1],v[-1]])
        ensemble_hierarchy.AppendToEnsemble(A,l)
    ensemble_hierarchy.EnsembleTransfer('Data')
    if np.shape(ensemble_hierarchy.Ensemble)!=(1,2,len(u[-1].dat.data),n):
        raise ValueError('failed')
    if ensemble_hierarchy.Type!='Data':
        raise ValueError('failed')


test_ensemble_type()



def test_ensemble_nxl():
    n=1
    L=1
    mesh=UnitSquareMesh(5,5)
    Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
    V=FunctionSpaceHierarchy(Mesh_Hierarchy,'DG',0)
    #
    if len(Mesh_Hierarchy)!=L+2:
        raise AssertionError('Mesh Hierarchy isnt of length L+2')
    #
    ensemble_hierarchy=EnsembleHierarchy()
    l=0
    for i in range(n):
        u=FunctionHierarchy(V)
        v=FunctionHierarchy(V)
        A=tuple([u[-1],v[-1]]) # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A,l)
    l=1
    for i in range(n):
        u=FunctionHierarchy(V)
        v=FunctionHierarchy(V)
        A=tuple([u[-1],v[-1]]) # tuple of same function space
        ensemble_hierarchy.AppendToEnsemble(A,l)
    # test nxl
    if hasattr(ensemble_hierarchy,'nxl')==0:
        raise AttributeError('failed. ensemble_hierarchy doesnt have nxl attribute')
    #
    if len(ensemble_hierarchy.nxl)!=2:
        raise ValueError('failed')
    if ensemble_hierarchy.nxl[0]!=(5*5):
        raise ValueError('failed')
    if ensemble_hierarchy.nxl[1]!=(10*10):
        raise ValueError('failed')

test_ensemble_nxl()





