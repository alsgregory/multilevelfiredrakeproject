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

from test_ensemble_case import *

ensemble_case_example = make_ensemble()

def test_discretization_timestepper(initial_condition_function,time_step_solve_function):
    mesh=UnitSquareMesh(10,10)
    L=1
    n_function_spaces=1
    vec=[0]
    families=['DG']
    degrees=[0]
    frame=setup(mesh,L,n_function_spaces,vec,families,degrees)
    frame.GenerateMeshHierarchy()
    frame.GenerateFunctionSpaceHierarchies()
    Mesh_Hierarchy=frame.Mesh_Hierarchy
    FunctionSpaceHierarchies=frame.FunctionSpaceHierarchies
    Courant=1.0
    lvlc=0
    lvlf=1
    nc=int(10*2**lvlc)
    nf=int(10*2**lvlf)
    hc=Timestep(nc,Courant).FindStableTimestep()
    hf=Timestep(nf,Courant).FindStableTimestep()
    sample = Discretization(lvlc,lvlf,hc,hf,initial_condition_function,time_step_solve_function,Mesh_Hierarchy,FunctionSpaceHierarchies)
    T=5
    sample.IC()
    sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc), we can now, then continue with state
    if abs(sample.solution.time-T)>1e-5:
        raise ValueError('failed')
    if len(sample.solution.state)!=2:
        raise ValueError('failed')


test_discretization_timestepper(initial_condition_function,time_step_solve_function)


def test_discretization_get_level(initial_condition_function,time_step_solve_function):
    mesh=UnitSquareMesh(10,10)
    L=1
    n_function_spaces=1
    vec=[0]
    families=['DG']
    degrees=[0]
    frame=setup(mesh,L,n_function_spaces,vec,families,degrees)
    frame.GenerateMeshHierarchy()
    frame.GenerateFunctionSpaceHierarchies()
    Mesh_Hierarchy=frame.Mesh_Hierarchy
    FunctionSpaceHierarchies=frame.FunctionSpaceHierarchies
    Courant=1.0
    lvlc=1
    lvlf=2
    nc=int(10*2**lvlc)
    nf=int(10*2**lvlf)
    hc=Timestep(nc,Courant).FindStableTimestep()
    hf=Timestep(nf,Courant).FindStableTimestep()
    sample = Discretization(lvlc,lvlf,hc,hf,initial_condition_function,time_step_solve_function,Mesh_Hierarchy,FunctionSpaceHierarchies)
    T=5
    sample.IC()
    sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc), we can now, then continue with state
    if get_level(sample.solution.state[0])[1]!=lvlc:
        raise ValueError('failed')
    if get_level(sample.solution.state[1])[1]!=lvlf:
        raise ValueError('failed')


test_discretization_get_level(initial_condition_function,time_step_solve_function)







