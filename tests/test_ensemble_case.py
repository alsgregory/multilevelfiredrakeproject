from __future__ import division # Get proper divison
import numpy as np
import random

from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level

from test_problem_functions import *

'''  Make sample ensembles for rest of test cases '''


def make_ensemble():
    mesh=UnitSquareMesh(4,4)
    L=2
    # Generate the Mesh Hierarchy and the FunctionSpaceHierarchies
    Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
    FunctionSpaceHierarchies=FunctionSpaceHierarchy(Mesh_Hierarchy,'CG',1)
    # Create Discretization Object
    ensemble_hierarchy=EnsembleHierarchy()
    level_to_prolong_to=3
    Ns=(64*2**(-np.linspace(0,2,3))).astype(int)
    for i in range(3):
        lvlc=i
        lvlf=i+1
        for j in range(Ns[i]):
            u_h=FunctionHierarchy(FunctionSpaceHierarchies)
            u=state(u_h[lvlc],u_h[lvlf])
            xc=np.random.normal(0.0,0.05,1)
            u.state[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
            u.state[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
            u.state[0].assign(Solve(lvlc,FunctionSpaceHierarchies,u.state[0]))
            u.state[1].assign(Solve(lvlf,FunctionSpaceHierarchies,u.state[1]))
            u_new=state(ProlongUpToFinestLevel(u.state[0],FunctionHierarchy(FunctionSpaceHierarchies)),ProlongUpToFinestLevel(u.state[1],FunctionHierarchy(FunctionSpaceHierarchies)))
            ensemble_hierarchy.AppendToEnsemble(u_new.state,i)
            assert get_level(u_new.state[1])[1]==len(Mesh_Hierarchy)-1
    return ensemble_hierarchy

if __name__ == "__main__":
    import os
    import pytest
    pytest.main(os.path.abspath(__file__))
