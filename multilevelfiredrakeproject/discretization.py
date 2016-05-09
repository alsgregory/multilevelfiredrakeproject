"""


Discretization and Implementation


Alastair Gregory 2016


"""


# NOTE UPDATE: NEED TO GET RID OF ANY PALMING LEVEL_COARSE AND LEVEL_FINE IN -> NEED TO USE THE GET_LEVEL FUNCTION THATS NOW IMPORTED INTO PACKAGES THAT CAN BE USED TO GET LEVEL AT ANY TIME


from packages import *

class ProblemSet():
    
    """Takes in 3 user defined :class:`Function' that defines the uncertainty quantification problem.
    
    	:arg initial_condition_function: an initialization function that outputs a :class:`state'.
    		Subarguments
    		:arg lvlc:
    		:arg :class:`MeshHierarchy':
    		:arg FunctionSpaceHierarchies: User-defined. 
        
        :arg time_step_solve_function: a timestep function that outputs a :class:`state'
        	Subarguments
        	:arg :class:`state':
        	:arg MeshHierarchy:
        	:arg FunctionSpaceHierarchies: User-defined. 
        
        :arg quantity_of_interest: a quantity of interest function that outputs a :class:`state.prepared_state'
        	Subarguments
        	:arg :class:`state':
        	:arg lvl_to_prolong_to: This is the level that :class:`state.prepared_state' will be on
        	:arg desired_family:
        	:arg desired_degree:
        	:arg :class:`MeshHierarchy':
        	:arg FunctionSpaceHierarchies: User-defined
        
        """
    
    def __init__(self , initial_condition_function , time_step_solve_function , quantity_of_interest):
        self.initial=initial_condition_function
        self.timestep=time_step_solve_function
        self.QoI=quantity_of_interest



class Discretization():
    
    """Class to manage the discretization of one sample of a :class:`state' using a user defined :class:`ProblemSet'
    
    	:arg lvlc:
    	:arg :class:`ProblemSet':
    	:arg FunctionSpaceHierarchies:
    	:arg Courant: (optional) Courant number, default 1.0
    """
    
    def __init__(self,lvlc,problemset,FunctionSpaceHierarchies,Courant=1.0):
        self.Mesh_Hierarchy=FunctionSpaceHierarchies[0]._mesh_hierarchy
        n00=(self.Mesh_Hierarchy[0].topology.num_cells()/2.0)**(1.0/float(self.Mesh_Hierarchy[0].topology.cell_dimension()))
        self.Courant=Courant
        self.nc=float(int(n00))*(2**(lvlc)); self.nf=2*self.nc
        self.hc=self.FindStableTimestep(self.nc,self.Courant)
        self.hf=0.5*self.hc
        self.lvlc=lvlc; self.lvlf=self.lvlc+1
        if (self.lvlf>len(self.Mesh_Hierarchy)-1) or (self.lvlc>len(self.Mesh_Hierarchy)-1):
            raise IndexError('Cannot discretize with level greater than the Mesh Hierarchy length')
        if (self.lvlc<0) or (self.lvlf<0):
            raise ValueError('Cannot have a level -1. Coarsest Level is indexed with 0.')
        self.initial_condition_function=problemset.initial; self.time_step_solve_function = problemset.timestep; self.QoI=problemset.QoI
        self.FunctionSpaceHierarchies=FunctionSpaceHierarchies
    
    def Timestepper(self,T): # thus can repeat this instance over and over again with increments of T
        
        """This discretizes the :class:`state' until time T
        
        	:arg T:
        """
        
        # check that initialization has been done
        if hasattr(self.solution,'state')==0:
            raise AttributeError('Discretization Object hasnt been initialized')
        #
        start_t=np.copy(self.solution.time)
        Nt=int(T/self.hc)
        for _ in range(Nt): # THIS HAS BEEN CHANGED!! TIMESTEP WAS ONE TOO MANY - INCONSISTENT!
            self.solution.time+=self.hc # update time
            self.solution=self.time_step_solve_function(self.solution,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies)
    
    def QuantityOfInterest(self,desired_family,desired_degree,lvl_to_prolong_to=None,index_of_state=0):
        
        """This computes the quantity of interest given a user defined QoI in ProblemSet
        
        	:arg desired_family: The family of the finite element of the desired quantity of interest :class:`FunctionSpace'
        	:arg desired_degree: The degree of the finite element of the desired quantity of interest :class`FunctionSpace'
        	:arg lvl_to_prolong_to: (optional) default, finest level
        	:arg index_of_state: (optional) default, 0. Given multiple :class:`Functions' in the :class:`state' this indicates which index of the list the quantity of interest is.
        """ 
        
        # own defined default
        if lvl_to_prolong_to==None:
            lvl_to_prolong_to=len(self.Mesh_Hierarchy)-1 # if not specified, go to last level
        prepared_state=self.QoI(self.solution,lvl_to_prolong_to,desired_family,desired_degree,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies,index_of_state)
        self.solution.prepared_state=prepared_state.state
    
    def IC(self):
        
        """Initializes a :class:`state' ready for discretization
        """
        
        self.solution=self.initial_condition_function(self.lvlc,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies)
        self.solution.lvlf=self.lvlf; self.solution.hf = self.hf
        self.solution.lvlc=self.lvlc; self.solution.hc = self.hc
        setattr(self.solution,'time',0.0)
    
    def FindStableTimestep(self,MeshPoints,Courant):
        """ Finds timstep satisfying stable courant number, given number of cells
        
        	:arg MeshPoints: Number of cells on that level mesh
        	:arg Courant: Courant Number
        
        """
        return (Courant / (MeshPoints))


