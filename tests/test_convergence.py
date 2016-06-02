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


def test_convergence_decerasing_nl():
    mesh=UnitSquareMesh(4,4)
    L=4
    # Generate the Mesh Hierarchy and the FunctionSpaceHierarchies
    Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
    FunctionSpaceHierarchies=FunctionSpaceHierarchy(Mesh_Hierarchy,'CG',1)
    # Create Discretization Object
    ensemble_hierarchy=EnsembleHierarchy()
    level_to_prolong_to=5
    eps=4e-2
    n=16
    converge=0
    i=0
    Nl=[]
    while converge==0:
        lvlc=i
        lvlf=i+1
        for j in range(n):
            u_h=FunctionHierarchy(FunctionSpaceHierarchies)
            u=state(u_h[lvlc],u_h[lvlf])
            xc=np.random.normal(0.0,0.05,1)
            u.state[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
            u.state[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
            u.state[0].assign(Solve(lvlc,FunctionSpaceHierarchies,u.state[0]))
            u.state[1].assign(Solve(lvlf,FunctionSpaceHierarchies,u.state[1]))
            u_new=state(ProlongUpToFinestLevel(u.state[0],FunctionHierarchy(FunctionSpaceHierarchies)),ProlongUpToFinestLevel(u.state[1],FunctionHierarchy(FunctionSpaceHierarchies)))
            ensemble_hierarchy.AppendToEnsemble(u_new.state,i)
            if get_level(u_new.state[1])[1]!=len(Mesh_Hierarchy)-1:
                raise ValueError('Prepared state does not have default level_to_prolong_to')
        Nl.append(n)
        SampleStatistics(ensemble_hierarchy)
        Ns=OptimalNl(ensemble_hierarchy,eps)
        for j in range(i+1):
            if Ns[j]>Nl[j]:
                dN=Ns[j]-Nl[j]
                for k in range(dN):
                    u_h=FunctionHierarchy(FunctionSpaceHierarchies)
                    u=state(u_h[lvlc],u_h[lvlf])
                    xc=np.random.normal(0.0,0.05,1)
                    u.state[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
                    u.state[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
                    u.state[0].assign(Solve(lvlc,FunctionSpaceHierarchies,u.state[0]))
                    u.state[1].assign(Solve(lvlf,FunctionSpaceHierarchies,u.state[1]))
                    u_new=state(ProlongUpToFinestLevel(u.state[0],FunctionHierarchy(FunctionSpaceHierarchies)),ProlongUpToFinestLevel(u.state[1],FunctionHierarchy(FunctionSpaceHierarchies)))
                    ensemble_hierarchy.AppendToEnsemble(u_new.state,i)
                Nl[j]=Ns[j]
        SampleStatistics(ensemble_hierarchy)
        converge=Convergence(ensemble_hierarchy,eps)
        i+=1
    if np.all(np.diff(np.array(Nl))<=0)==0:
        raise AssertionError('failed. hasnt converged with nl decreasing')
    return ensemble_hierarchy

test_convergence_decerasing_nl()


