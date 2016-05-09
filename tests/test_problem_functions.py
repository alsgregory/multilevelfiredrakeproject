"""

Case Study User Defined Functions -> Heat Equation


Alastair Gregory 2016

"""

from __future__ import division # Get proper divison
import numpy as np
import random
from scipy import stats
from scipy.stats import norm
from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level

def initial_condition_function(level_coarse,level_fine,Mesh,FunctionSpaceHierarchies): # CoarseLevelIndex starts at 0 for L=0 # Random x_spread and y_spread!
    VcgH = FunctionSpaceHierarchies[0]
    #VuH = VectorFunctionSpaceHierarchy(Mesh,"DG",1)
    u=FunctionHierarchy(VcgH)
    # now make individual 'sub' hierarchies
    uh=[u[level_coarse],u[level_fine]]
    # Initial Condition
    xc=np.random.normal(0.0,0.25,1)
    uh[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    uh[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    # generate a state object here
    State=state(uh[0],uh[1])
    return State

def f(u,v,lvl,h,BM):
    Vcg = v
    v = TestFunction(Vcg[lvl])
    u_ = TrialFunction(Vcg[lvl])
    # get information about the level, timestep and dx
    dt = h
    a=(dt*dot(grad(v),grad(u_)) + v*u_)*dx # remember to multiply Laplacian by mu and dt
    L=(u+BM)*v*dx # RHS is function on RHS of Helmholtz
    u_new=FunctionHierarchy(Vcg)[lvl]
    u_problem = LinearVariationalProblem(a,L,u_new)
    solve(a==L,u_new)
    return u_new

def GenerateRandomGaussianField(FunctionSpace_Hierarchy,fine_level,df,xi):
    dWC=FunctionHierarchy(FunctionSpace_Hierarchy)[fine_level-1].assign(0)
    for i in range(2):
        dW=FunctionHierarchy(FunctionSpace_Hierarchy)
        dW[fine_level].assign(0)
        dW[fine_level-1].assign(0)
        dwf=np.sqrt(df)*xi.dat.data*np.random.normal(0,1,np.shape(dW[fine_level].dat.data))
        dW[fine_level].dat.data[:]=dwf
        Comparison=FunctionHierarchy(FunctionSpace_Hierarchy)[fine_level-1].assign(0)
        inject(dW[fine_level],Comparison)
        dWC.dat.data[:] += Comparison.dat.data
        if i==0:
            dW2=FunctionHierarchy(FunctionSpace_Hierarchy)
            dW2[fine_level].assign(0)
            dW2[fine_level].dat.data[:]=np.copy(dW[fine_level].dat.data)
        else:
            dW[fine_level].dat.data[:]=dW[fine_level].dat.data
    return([dWC,[dW2[fine_level],dW[fine_level]]])

def time_step_solve_function(State,Mesh,FunctionSpaceHierarchies):
    """
    Hierarchy sets are from last time step. They are the functions as follows:
    dq1,qh,q1,psi0,psi1,q0. Function spaces can be gained from the original heirarchy set ups.
    """
    # Get Function spaces and state
    uh=State.state
    Vcg=FunctionSpaceHierarchies[0]
    lvlc = State.levels[0]
    lvlf = State.levels[1]
    # u derivative
    dtc = State.hc
    dtf = State.hf
    xi=Constant(0.25)
    BM=GenerateRandomGaussianField(Vcg,lvlf,dtf,xi)
    for j in range(2):
        #Predictor stage
        #u=uh[1]
        u_new=f(uh[1],Vcg,lvlf,dtf,BM[1][j])
        uh[1].assign(u_new)
    """
    Coarse Solve
    """
    # Predictor stage
    #u=uh[0]
    u_new=f(uh[0],Vcg,lvlc,dtc,BM[0])
    uh[0].assign(u_new)
    return State


def QoI(solution,lvl_to_prolong_to,desired_family,desired_degree,Mesh_Hierarchy,FunctionSpaceHierarchies,index_of_state=0):
    # create comparison function hierarchies
    FHc=FunctionHierarchy(FunctionSpaceHierarchy(Mesh_Hierarchy,desired_family,desired_degree))[lvl_to_prolong_to]
    FHf=FunctionHierarchy(FunctionSpaceHierarchy(Mesh_Hierarchy,desired_family,desired_degree))[lvl_to_prolong_to]
    lvlc=solution.lvlc
    lvlf=solution.lvlf
    # prolong these two
    deg=solution.state[0].function_space().ufl_element().degree()
    family=solution.state[0].function_space().ufl_element().family()
    CompCoarse=FunctionHierarchy(FunctionSpaceHierarchies[0])
    CompFine=FunctionHierarchy(FunctionSpaceHierarchies[0])
    FHc.project(ProlongInject().ProlongUpToAnyLevel(lvlc,lvl_to_prolong_to,solution.state[0],CompCoarse))
    FHf.project(ProlongInject().ProlongUpToAnyLevel(lvlf,lvl_to_prolong_to,solution.state[1],CompFine))
    solution_copy=solution
    new_tuple=tuple([FHc,FHf])
    solution_copy.state=new_tuple
    return solution_copy


