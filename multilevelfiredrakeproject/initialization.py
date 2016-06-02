

from packages import *


def GenerateMeshHierarchy(coarsest_mesh,L):
    
    """ Set ups the :class:`MeshHierarchy' needed for a mlmc implementation for a given L. This will have L+2 possible levels of resolution. l=-1,...,L.
    
    	:param coarsest_mesh: The coarsest mesh, level -1 of the hierarchy of levels.
    	:type coarsest_mesh: :class:`Mesh'
    	
    	:param L: L
    	:type L: int
    
    """
    
    hierarchy=MeshHierarchy(coarsest_mesh,L+1)
    Mesh_Hierarchy=hierarchy
    return Mesh_Hierarchy
    


