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
    mesh=UnitSquareMesh(2,2)
    L=2
    n_function_spaces=1
    vec=[0]
    families=['DG']
    degrees=[1]
    frame=setup(mesh,L,n_function_spaces,vec,families,degrees)
    # Generate the Mesh Hierarchy and the FunctionSpaceHierarchies
    frame.GenerateMeshHierarchy()
    frame.GenerateFunctionSpaceHierarchies()
    Mesh_Hierarchy=frame.Mesh_Hierarchy
    FunctionSpaceHierarchies=frame.FunctionSpaceHierarchies
    # Create Discretization Object
    T=0.25
    Courant=0.25
    ensemble_hierarchy=EnsembleHierarchy()
    deg=1
    fam='DG'
    level_to_prolong_to=3
    Ns=(32*2**(-np.linspace(0,2,3))).astype(int)
    for i in range(3):
        lvlc=i
        lvlf=i+1
        nc=int(10*2**lvlc)
        nf=int(10*2**lvlf)
        hc=Timestep(nc,Courant).FindStableTimestep()
        hf=Timestep(nf,Courant).FindStableTimestep()
        for j in range(Ns[i]):
            sample = Discretization(lvlc,lvlf,hc,hf,initial_condition_function,time_step_solve_function,Mesh_Hierarchy,FunctionSpaceHierarchies)
            sample.IC()
            sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc),
            sample.QuantityOfInterest(QoI,level_to_prolong_to,fam,deg,Mesh_Hierarchy)
            ensemble_hierarchy.AppendToEnsemble(sample.solution,i) # try and do this automatically depending on the level (i.e. know where to append into) - given in sample and solution
    return ensemble_hierarchy



