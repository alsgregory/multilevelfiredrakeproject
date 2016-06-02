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




mesh=UnitSquareMesh(5,5)

L=2

Mesh_Hierarchy=GenerateMeshHierarchy(mesh,L)
FunctionSpaceHierarchies=FunctionSpaceHierarchy(Mesh_Hierarchy,'CG',1)

# Create Discretization Object


def solve1(lvl,FunctionSpaceHierarchies,u):
    v=TestFunction(FunctionSpaceHierarchies[lvl])
    u_ = TrialFunction(FunctionSpaceHierarchies[lvl])
    a=(dot(grad(v),grad(u_)) + v*u_)*dx # remember to multiply Laplacian by mu and dt
    L=(u)*v*dx # RHS is function on RHS of Helmholtz
    u_new=FunctionHierarchy(FunctionSpaceHierarchies)[lvl]
    solve(a==L,u_new)
    return u_new



ensemble_hierarchy=EnsembleHierarchy()


lvlc=0
lvlf=1


for i in range(64):
    u_h=FunctionHierarchy(FunctionSpaceHierarchies)
    u=state(u_h[lvlc],u_h[lvlf])
    # initial condition
    xc=np.random.normal(0.0,0.05,1)
    u.state[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    u.state[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    # solve
    u.state[0].assign(solve1(lvlc,FunctionSpaceHierarchies,u.state[0]))
    u.state[1].assign(solve1(lvlf,FunctionSpaceHierarchies,u.state[1]))
    u_new=state(ProlongUpToFinestLevel(u.state[0],FunctionHierarchy(FunctionSpaceHierarchies)),ProlongUpToFinestLevel(u.state[1],FunctionHierarchy(FunctionSpaceHierarchies)))
    ensemble_hierarchy.AppendToEnsemble(u_new.state,0)


lvlc=1
lvlf=2


for i in range(32):
    u_h=FunctionHierarchy(FunctionSpaceHierarchies)
    u=state(u_h[lvlc],u_h[lvlf])
    # initial condition
    xc=np.random.normal(0.0,0.05,1)
    u.state[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    u.state[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    # solve
    u.state[0].assign(solve1(lvlc,FunctionSpaceHierarchies,u.state[0]))
    u.state[1].assign(solve1(lvlf,FunctionSpaceHierarchies,u.state[1]))
    u_new=state(ProlongUpToFinestLevel(u.state[0],FunctionHierarchy(FunctionSpaceHierarchies)),ProlongUpToFinestLevel(u.state[1],FunctionHierarchy(FunctionSpaceHierarchies)))
    ensemble_hierarchy.AppendToEnsemble(u_new.state,1)


lvlc=2
lvlf=3


for i in range(16):
    u_h=FunctionHierarchy(FunctionSpaceHierarchies)
    u=state(u_h[lvlc],u_h[lvlf])
    # initial condition
    xc=np.random.normal(0.0,0.05,1)
    u.state[0].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    u.state[1].interpolate(Expression("exp(-(pow((x[0]-0.5+x_center),2)/0.25)-(pow(x[1]-0.5,2)/0.25))",x_center=xc))
    # solve
    u.state[0].assign(solve1(lvlc,FunctionSpaceHierarchies,u.state[0]))
    u.state[1].assign(solve1(lvlf,FunctionSpaceHierarchies,u.state[1]))
    u_new=state(ProlongUpToFinestLevel(u.state[0],FunctionHierarchy(FunctionSpaceHierarchies)),ProlongUpToFinestLevel(u.state[1],FunctionHierarchy(FunctionSpaceHierarchies)))
    ensemble_hierarchy.AppendToEnsemble(u_new.state,2)













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

OptimalNl(ensemble_hierarchy,1e-3)
Convergence(ensemble_hierarchy,1e-3)

# Change function space via Forecast or something, and don't update SampleStats.

ensemble_hierarchy.EnsembleTransfer('Function')

SampleStatistics(ensemble_hierarchy)

# Then see if Bounds Have Error

OptimalNl(ensemble_hierarchy,1e-3)
Convergence(ensemble_hierarchy,1e-3)

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



