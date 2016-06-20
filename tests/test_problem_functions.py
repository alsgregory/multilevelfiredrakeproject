"""

Case Study User Defined Functions -> Heat Equation


Alastair Gregory 2016

"""

from __future__ import division # Get proper divison
import numpy as np
import random

from firedrake import *
parameters["reorder_meshes"] = False
from multilevelfiredrakeproject import *
from firedrake.mg.utils import get_level


def Solve(lvl,FunctionSpaceHierarchies,u):
    v=TestFunction(FunctionSpaceHierarchies[lvl])
    u_ = TrialFunction(FunctionSpaceHierarchies[lvl])
    a=(dot(grad(v),grad(u_)) + v*u_)*dx # remember to multiply Laplacian by mu and dt
    L=(u)*v*dx # RHS is function on RHS of Helmholtz
    u_new=FunctionHierarchy(FunctionSpaceHierarchies)[lvl]
    solve(a==L,u_new)
    return u_new





