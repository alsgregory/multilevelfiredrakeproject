

from packages import *

class setup():
    def __init__(self,coarsest_mesh,L,FunctionSpaces):
        """ FunctionSpaces is a list of all FS used. The mesh corresponds to coarsest_mesh. This is used to extract families and degrees / vector info """
        self.coarsest_mesh=coarsest_mesh; self.L=L
        self.FunctionSpaces=FunctionSpaces
        if type(self.FunctionSpaces)==FunctionSpace: # if just one FS
            self.FunctionSpaces=[self.FunctionSpaces]
        self.n_function_spaces=len(self.FunctionSpaces)
        # now generate families and degrees from FunctionSpaces list.
        self.vec=[]
        self.families=[]
        self.degrees=[]
        for i in range(self.n_function_spaces):
            v=float(self.FunctionSpaces[i].ufl_element().value_size())-1.0
            d=self.FunctionSpaces[i].ufl_element().degree()
            f=self.FunctionSpaces[i].ufl_element().family()
            if f=='Lagrange':
                ff='CG'
            if f=='Discontinuous Lagrange':
                ff='DG'
            self.vec.append(v)
            self.families.append(ff)
            self.degrees.append(d)
    def GenerateMeshHierarchy(self): # mesh hierarchy will be L+1, so do finest level L here
        """ Creates a mesh hierarchy of all the levels, from mesh on coarsest level
        Outputs:
        - hierarchy, mesh hierarchy
        For a list of available firedrake meshes, go to http://firedrakeproject.org/variational-problems.html#constructing-meshes
        """
        hierarchy=MeshHierarchy(self.coarsest_mesh,self.L+1)
        self.Mesh_Hierarchy=hierarchy
    def GenerateFunctionSpaceHierarchies(self):
        """ Creates a list of all function space hierarchies. vec is a list, size n_function_spaces, of 1,0 if a vector, familes and degrees self explanatory """
        if hasattr(self,'Mesh_Hierarchy')==0:
            raise AttributeError('setup doesnt have attribute Mesh_Hierarchy')
        l=[]
        for i in range(self.n_function_spaces):
            if self.vec[i]==1:
                l.append(VectorFunctionSpaceHierarchy(self.Mesh_Hierarchy,self.families[i],self.degrees[i]))
            if self.vec[i]==0:
                l.append(FunctionSpaceHierarchy(self.Mesh_Hierarchy,self.families[i],self.degrees[i]))
        self.FunctionSpaceHierarchies=l


