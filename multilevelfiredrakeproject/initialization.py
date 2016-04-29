"""


Initialization


Alastair Gregory 2016

"""


from packages import *

class setup():
    def __init__(self,coarsest_mesh,L,n_function_spaces,vec,families,degrees):
        self.coarsest_mesh=coarsest_mesh; self.L=L
        self.n_function_spaces=n_function_spaces
        self.vec=vec
        self.families=families
        self.degrees=degrees
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


