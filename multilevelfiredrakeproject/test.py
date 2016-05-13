"""


test - This tests out all types, structures and modules housed by the software - multilevelfiredrakeproject

- with the Heat Equation -> gives an ensemble forecast gained by MLMC ensemble forecasting framework

- tests out types -> EnsembleHierarchy, EnsembleForecast function, State


"""


from __future__ import division # Get proper divison

import itertools
import numpy as np
import math
import random 
from scipy import integrate as int1
from scipy.stats import gaussian_kde
from scipy import stats
import pickle
from scipy.stats import norm
import time

from firedrake import *
parameters["reorder_meshes"] = False

from firedrake.mg.utils import get_level
from multilevelfiredrakeproject import *


def initial_condition_function(level_coarse,Mesh,FunctionSpaceHierarchies): # CoarseLevelIndex starts at 0 for L=0 # Random x_spread and y_spread!
    VcgH = FunctionSpaceHierarchies[0]
    #VuH = VectorFunctionSpaceHierarchy(Mesh,"DG",1)
    u=FunctionHierarchy(VcgH)
    level_fine=level_coarse+1
    # now make individual 'sub' hierarchies
    uh=[u[level_coarse],u[level_fine]]
    # Initial Condition
    xc=np.random.normal(0.0,0.05,1)
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
    perp=lambda u: as_vector((u,u))
    a=(dt*dot(grad(v),grad(u_)) + v*u_)*dx - Constant(0.5)*dt*dot(perp(v),grad(u_))*dx # remember to multiply Laplacian by mu and dt
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
    xi=Constant(0.005)
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
    FHc.project(ProlongInject().ProlongUpToAnyLevel(lvl_to_prolong_to,solution.state[0],CompCoarse))
    FHf.project(ProlongInject().ProlongUpToAnyLevel(lvl_to_prolong_to,solution.state[1],CompFine))
    solution_copy=solution
    new_tuple=tuple([FHc,FHf])
    solution_copy.state=new_tuple
    return solution_copy



problemset=ProblemSet(initial_condition_function,time_step_solve_function,QoI)


mesh=UnitSquareMesh(2,2)

L=2

FunctionSpaces=FunctionSpace(mesh,'CG',1)

frame=setup(mesh,L,FunctionSpaces)

# Generate the Mesh Hierarchy and the FunctionSpaceHierarchies

frame.GenerateMeshHierarchy()
frame.GenerateFunctionSpaceHierarchies()

Mesh_Hierarchy=frame.Mesh_Hierarchy
FunctionSpaceHierarchies=frame.FunctionSpaceHierarchies

# Create Discretization Object

Courant=0.0625

lvlc=0
lvlf=1
nc=int(10*2**lvlc)
nf=int(10*2**lvlf)



T=0.0625

ensemble_hierarchy=EnsembleHierarchy()

deg=1
fam='CG'
level_to_prolong_to=3



for i in range(512):
    sample = Discretization(lvlc,problemset,FunctionSpaceHierarchies,Courant)
    sample.IC()
    sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc), we can now, then continue with state
    sample.QuantityOfInterest(fam,deg)
    ensemble_hierarchy.AppendToEnsemble(sample.solution,0)
    print sample.hf
    print sample.hc
    print sample.lvlc
    print sample.lvlf



lvlc=1
lvlf=2
nc=int(10*2**lvlc)
nf=int(10*2**lvlf)



for i in range(256):
    sample = Discretization(lvlc,problemset,FunctionSpaceHierarchies,Courant)
    sample.IC()
    sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc), we can now, then continue with state
    sample.QuantityOfInterest(fam,deg)
    ensemble_hierarchy.AppendToEnsemble(sample.solution,1) # try and do this automatically depending on the level (i.e. know where to append into) - given in sample and solution


lvlc=2
lvlf=3
nc=int(10*2**lvlc)
nf=int(10*2**lvlf)



for i in range(128):
    sample = Discretization(lvlc,problemset,FunctionSpaceHierarchies,Courant)
    sample.IC()
    sample.Timestepper(T) # if we wanted to do something to state (importance sampling etc), we can now, then continue with state
    sample.QuantityOfInterest(fam,deg)
    ensemble_hierarchy.AppendToEnsemble(sample.solution,2)







# check

print ensemble_hierarchy.Type

SampleStatistics(ensemble_hierarchy)

print 'mean: ', norm(ensemble_hierarchy.Mean[0]), norm(ensemble_hierarchy.Mean[1]), norm(ensemble_hierarchy.Mean[2])
print 'variance: ', norm(ensemble_hierarchy.Variance[0]), norm(ensemble_hierarchy.Variance[1]), norm(ensemble_hierarchy.Variance[2])


MLMCMean=ensemble_hierarchy.MultilevelExpectation

# CHANGE to check stats work in other type

ensemble_hierarchy.EnsembleTransfer('Data')

print ensemble_hierarchy.Type

SampleStatistics(ensemble_hierarchy)

print 'mean: ', np.linalg.norm(ensemble_hierarchy.Mean[0]),np.linalg.norm(ensemble_hierarchy.Mean[1]), np.linalg.norm(ensemble_hierarchy.Mean[2])
print 'variance: ', np.linalg.norm(ensemble_hierarchy.Variance[0]),np.linalg.norm(ensemble_hierarchy.Variance[1]), np.linalg.norm(ensemble_hierarchy.Variance[2])

MLMCMean=ensemble_hierarchy.MultilevelExpectation





""" ********** """

# find error bounding


Bounds=ErrorBounding(ensemble_hierarchy,1e-3)
Bounds.OptimalNl()
Bounds.Convergence()

# Change function space via Forecast or something, and don't update SampleStats.

ensemble_hierarchy.EnsembleTransfer('Function')

SampleStatistics(ensemble_hierarchy)

# Then see if Bounds Have Error

Bounds=ErrorBounding(ensemble_hierarchy,1e-3)
Bounds.OptimalNl()
Bounds.Convergence()

""" ********** """

# CREATE ENSEMBLE FORECAST - the ensemble hierarchy shouldn't change



Weights=ensemble_hierarchy.Weights
N=256
loc=False
Sigma=100

Forecast=EnsembleForecast(ensemble_hierarchy,ensemble_hierarchy._EnsembleHierarchy__OriginalFunctionSpaces[-1][1],N,loc,Sigma)

Forecast.EnsembleTransfer('Data')

# converted back to data

#############################################

ensemble_hierarchy.EnsembleTransfer('Function')

C=Forecast.Cov

# State what delta used to shrink covariance matrix

delta=Forecast.Delta


##############################################


Forecast.EnsembleTransfer('Function')

FileMember = File("ensemble_member_q.pvd")
Forecast.Forecast[3].rename('function')
FileMember << Forecast.Forecast[5]

ensemble_hierarchy.EnsembleTransfer('Function')
MLMCMean=ensemble_hierarchy.MultilevelExpectation
# Multilevel Mean

FileMember = File("mean_q.pvd")
MLMCMean.rename('function')
FileMember << MLMCMean


# List what attrbiutes an object has. e.g. ensemble_hierarchy

dir(ensemble_hierarchy)

# coverage

c=0
for i in range(N):
    if Forecast.Forecast[i].at([1,0])>0.55 and Forecast.Forecast[i].at([0,1])<0.55:
        c+=1
        print c


ensemble_hierarchy.EnsembleTransfer('Function')

c_standard=0
for i in range(len(ensemble_hierarchy.Ensemble[-1])):
    if ensemble_hierarchy.Ensemble[-1][i][1].at([1,0])>0.55 and ensemble_hierarchy.Ensemble[-1][i][1].at([0,1])<0.55:
        c_standard+=1
        print c_standard


p_c=c/N

p_c_standard=c_standard/len(ensemble_hierarchy.Ensemble[-1])



