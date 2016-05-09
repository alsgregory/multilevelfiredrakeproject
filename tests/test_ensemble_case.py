from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level

from test_problem_functions import *

'''  Make sample ensembles for rest of test cases '''


def make_ensemble():
    mesh=UnitSquareMesh(4,4)
    L=2
    FunctionSpaces=FunctionSpace(mesh,'CG',1)
    frame=setup(mesh,L,FunctionSpaces)
    # Generate the Mesh Hierarchy and the FunctionSpaceHierarchies
    frame.GenerateMeshHierarchy()
    frame.GenerateFunctionSpaceHierarchies()
    Mesh_Hierarchy=frame.Mesh_Hierarchy
    FunctionSpaceHierarchies=frame.FunctionSpaceHierarchies
    # Create Discretization Object
    T=0.25
    Courant=0.25
    ensemble_hierarchy=EnsembleHierarchy()
    problemset=ProblemSet(initial_condition_function,time_step_solve_function,QoI)
    deg=1
    fam='CG'
    level_to_prolong_to=3
    Ns=(64*2**(-np.linspace(0,2,3))).astype(int)
    for i in range(3):
        lvlc=i
        lvlf=i+1
        for j in range(Ns[i]):
            sample = Discretization(lvlc,problemset,FunctionSpaceHierarchies,Courant)
            sample.IC()
            sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc),
            sample.QuantityOfInterest(fam,deg)
            if get_level(sample.solution.prepared_state[1])[1]!=len(Mesh_Hierarchy)-1:
                raise ValueError('Prepared state does not have default level_to_prolong_to')
            ensemble_hierarchy.AppendToEnsemble(sample.solution,i) # try and do this automatically depending on the level (i.e. know where to append into) - given in sample and solution
    return ensemble_hierarchy



