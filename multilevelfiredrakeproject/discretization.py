

from packages import *

class ProblemSet():
    
    """Takes in 3 user defined :attr:`function` that defines the uncertainty quantification problem.
    
    	:param initial_condition_function: An initialization function that only outputs a :class:`state`.
    	
    	:type initial_condition_function: function, args: lvlc,:class:`MeshHierarchy`,FunctionSpaceHierarchies
    	
    		:param FunctionSpaceHierarchies: A list of multiple :class:`FunctionSpaceHierarchy`.
    		:type FunctionSpaceHierarchies: list
    	
    	:param time_step_solve_function: A timestep function that only outputs a :class:`state`.
    	:type time_step_solve_function: function, args: :class:`state`,:class:`MeshHierarchy`,FunctionSpaceHierarchies
    	
    	:param quantity_of_interest: A quantity of interest function that only outputs a :attr:`prepared_state'
    	:type quantity_of_interest: function(:class:`state`,lvl_to_prolong_to,desired_family,desired_degree,:class:`MeshHierarchy`,FunctionSpaceHierarchies)
    		
    		:param lvl_to_prolong_to: Level to prolong functions from different levels to in quantity of interest. e.g. typically the finest level.
    		:type lvl_to_prolong_to: int
    		
    		:param desired_family: Family of Finite Element for the desired :class:`FunctionSpace` for the quantity of interest :class:`Function` to exist on.
    		:type desired_family: string
    		
    		:param desired_degree: Degree of Finite Element for the desired :class:`FunctionSpace` for the quantity of interest :class:`Function` to exist on.
    		:type desired_degree: string
        
        """
    
    def __init__(self , initial_condition_function , time_step_solve_function , quantity_of_interest):
        self.initial=initial_condition_function
        self.timestep=time_step_solve_function
        self.QoI=quantity_of_interest



class Discretization():
    
    """Class to manage the discretization of one sample of a :class:`state' using a user defined :class:`ProblemSet'.
    
    	:arg lvlc:
    	:arg :class:`ProblemSet':
    	:arg FunctionSpaceHierarchies:
    	:arg Courant: (optional) default=1.0. Courant number.
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
        
        """This discretizes the :class:`state' until time T.
        
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
            # check for type
            if hasattr(self.solution,'state')!=1:
                raise TypeError('Output from time_step_solve_function is not of state type')
    
    def QuantityOfInterest(self,desired_family,desired_degree,lvl_to_prolong_to=None,index_of_state=0):
        
        """This computes the quantity of interest given a user defined QoI in ProblemSet.
        
        	:arg desired_family: The family of the finite element of the desired quantity of interest :class:`FunctionSpace'
        	:arg desired_degree: The degree of the finite element of the desired quantity of interest :class`FunctionSpace'
        	:arg lvl_to_prolong_to: (optional) default=finest level
        	:arg index_of_state: (optional) default=0. Given multiple :class:`Functions' in the :class:`state' this indicates which index of the list the quantity of interest is.
        	
        """ 
        
        # own defined default
        if lvl_to_prolong_to==None:
            lvl_to_prolong_to=len(self.Mesh_Hierarchy)-1 # if not specified, go to last level
        prepared_state=self.QoI(self.solution,lvl_to_prolong_to,desired_family,desired_degree,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies,index_of_state)
        self.solution.prepared_state=prepared_state.state
    
    def IC(self):
        
        """Initializes a :class:`state' ready for discretization.
        """
        
        self.solution=self.initial_condition_function(self.lvlc,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies)
        # Test the output
        if hasattr(self.solution,'state')!=1:
            raise TypeError('Output from initial_condition_function is not of state type')
        self.solution.lvlf=self.lvlf; self.solution.hf = self.hf
        self.solution.lvlc=self.lvlc; self.solution.hc = self.hc
        setattr(self.solution,'time',0.0)
    
    def FindStableTimestep(self,MeshPoints,Courant):
        """ Finds timstep satisfying stable courant number, given number of cells.
        
        	:param MeshPoints: Number of cells on that level mesh.
        	:type MeshPoints: int
        	
        	:param Courant: Courant Number.
        	:type Courant: float
        
        """
        return (Courant / (MeshPoints))


