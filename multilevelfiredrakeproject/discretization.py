"""


Discretization and Implementation


Alastair Gregory 2016


"""


# NOTE UPDATE: NEED TO GET RID OF ANY PALMING LEVEL_COARSE AND LEVEL_FINE IN -> NEED TO USE THE GET_LEVEL FUNCTION THATS NOW IMPORTED INTO PACKAGES THAT CAN BE USED TO GET LEVEL AT ANY TIME


from packages import *

class Discretization():
    def __init__(self,lvlc,lvlf,hc,hf,initial_condition_function,time_step_solve_function,Mesh_Hierarchy,FunctionSpaceHierarchies):
        self.lvlc=lvlc; self.lvlf=lvlf
        if (lvlf>len(Mesh_Hierarchy)-1) or (lvlc>len(Mesh_Hierarchy)-1):
            raise IndexError('Cannot discretize with level greater than the Mesh Hierarchy length')
        if (lvlc<0) or (lvlf<0):
            raise ValueError('Cannot have a level -1. Coarsest Level is indexed with 0.')
        self.hc=hc; self.hf=hf
        self.initial_condition_function=initial_condition_function; self.time_step_solve_function = time_step_solve_function
        self.Mesh_Hierarchy=Mesh_Hierarchy
        self.FunctionSpaceHierarchies=FunctionSpaceHierarchies
    def Timestepper(self,T): # thus can repeat this instance over and over again with increments of T
        """ Finds the coarse and fine functions after solving the function @time_step_solve_function one timestep (dtc and dtf) - at the same time so they have same forcing
        Inputs:
        - T, final time
        - hc, accuracy parameter (time-step) on finer level
        - hf, accuracy parameter (time-step) on coarse level
        - time_step_solve_function, function for each time step of coarse, with fine . Args: Tuple of 2 x tuple_solution, tuple of neccersary functions (Initial_Inputs), hc, hf (these accuracy parameters will have to be included to define levels), mesh hierarchy
        - Initial_Inputs, tuple of 2 x tuple containing the initial inputs into time_step_solve_function, can be at time 0, or last assimilation step
        Outputs:
        - FinalOuputs, tuple containing Initial_Inputs having gone through suitbale amount of timesteps
        """
        # check that initialization has been done
        if hasattr(self.solution,'state')==0:
            raise AttributeError('Discretization Object hasnt been initialized')
        #
        start_t=np.copy(self.solution.time)
        Nt=int(T/self.hc)
        for _ in range(Nt): # THIS HAS BEEN CHANGED!! TIMESTEP WAS ONE TOO MANY - INCONSISTENT!
            self.solution.time+=self.hc # update time
            self.solution=self.time_step_solve_function(self.solution,self.lvlc,self.lvlf,self.hc,self.hf,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies)
    def QuantityOfInterest(self,QoI,lvl_to_prolong_to,desired_family,desired_degree,index_of_state=0):
        prepared_state=QoI(self.solution,lvl_to_prolong_to,desired_family,desired_degree,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies,index_of_state)
        self.solution.prepared_state=prepared_state.state
    def IC(self):
        """ Initializes the state
        Inputs:
        - level_coarse, finer level
        - level_fine, coarse level
        - initial_condition, a function taking level_coarse and level_fine, initializing Initial_Intputs, tuple of 2 x tuple of Initial_Inputs for Timestepper
        Outputs:
        - Initial_Inputs
        """
        self.solution=self.initial_condition_function(self.lvlc,self.lvlf,self.Mesh_Hierarchy,self.FunctionSpaceHierarchies)
        self.solution.lvlf=self.lvlf; self.solution.hf = self.hf
        self.solution.lvlc=self.lvlc; self.solution.hc = self.hc
        setattr(self.solution,'time',0.0)


class Timestep():
    def __init__(self,MeshPoints,Courant):
        self.MeshPoints=MeshPoints; self.Courant=Courant
    def FindStableTimestep(self):
        """ Finds timstep satisfying stable courant number, given mesh points"""
        return (self.Courant / (self.MeshPoints))

