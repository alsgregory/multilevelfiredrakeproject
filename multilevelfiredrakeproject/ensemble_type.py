

from packages import *



class EnsembleHierarchy():
    """This class creates and builds a hierarchy of ensembles as well as containing properties of them, such as statistics.
    
    	:param M: Refinement factor.
    	:type M: int
    	
    """
    
    def __init__(self,M=2):
        self.Type="Function" #: This is either 'Data' or 'Function' depending on what the entries into :attr:`EnsembleHierarchy.Ensemble` are.
        self.Ensemble=[] #: This is the actual hierarchy of ensembles
        self.hl=[] #: A list of the finer timesteps used on each level of the hierarchy
        self.L=0 #: Length of the hierarchy
        self.M=M #: Refinement factor
        # Preallocate attributes
        self.MultilevelExpectation=None #: Either 'Data' or 'Function' of the MLMC estimator of the expectation of the ensembles.
        self.Mean=[] #: A list of the sample means of each of the different levels
        self.Variance=[] #: A list of the sample variances of each of the different levels
        #
        self.__OriginalFunctionSpaces=[]
        self.Weights=[] #: Same size as :attr:`EnsembleHierarchy.Ensemble` giving weights of the indidivual ensemble members. Default is evenly weighted.
    
    
    def __SameFunctionSpaceErrorCheck(self):
        
        Comparison=len(self.Ensemble[0][0][0].dat.data)
        for i in range(len(self.Ensemble)):
            for j in range(len(self.Ensemble[i])):
                for k in range(len(self.Ensemble[i][j])):
                    if len(self.Ensemble[i][j][k].dat.data)!=Comparison:
                        raise ValueError('All Ensemble Members Do Not Have Same Function Spaces')
    
    
    def AppendToEnsemble(self,tuple_to_append,Level_to_append_to):
        
        if hasattr(tuple_to_append,'prepared_state')==0:
            raise AttributeError('State hasnt been prepared for appedning into ensemble hierarchy')
        else:
            if tuple_to_append.prepared_state==None:
                raise ValueError('State hasnt been prepared for appending into ensemble hierarchy. Use quantity of interest function. See attribute discretization.ProblemSet.QoI')
        if self.L+1<=Level_to_append_to:
            raise ValueError('Cant append to level more than 1 above EnsembleHierarchy.L')
        #try:
        if Level_to_append_to<len(self.Ensemble):
            self.Ensemble[Level_to_append_to].append(tuple_to_append.prepared_state)
            self.__SameFunctionSpaceErrorCheck() # check for all same function space
            self.Weights[Level_to_append_to]=tuple([np.ones(len(self.Ensemble[Level_to_append_to]))*(1/float(len(self.Ensemble[Level_to_append_to]))),np.ones(len(self.Ensemble[Level_to_append_to]))*(1/float(len(self.Ensemble[Level_to_append_to])))])
        #except IndexError:
        if Level_to_append_to>=len(self.Ensemble):
            self.Ensemble.append([])
            self.Weights.append([])
            self.Ensemble[Level_to_append_to].append(tuple_to_append.prepared_state)
            self.L=self.L+1
            self.hl.append(tuple_to_append.hf)
            self.__OriginalFunctionSpaces.append(tuple([tuple_to_append.prepared_state[0].function_space(),tuple_to_append.prepared_state[1].function_space()]))
            self.__SameFunctionSpaceErrorCheck() # check for all same function space
            self.Weights[Level_to_append_to]=tuple([np.ones(len(self.Ensemble[Level_to_append_to]))*(1/float(len(self.Ensemble[Level_to_append_to]))),np.ones(len(self.Ensemble[Level_to_append_to]))*(1/float(len(self.Ensemble[Level_to_append_to])))])
    
    
    def EnsembleTransfer(self,To="Data"):
        
        if self.Type!=To: # only change type if not already that
            if To=="Data":
                self.__SameFunctionSpaceErrorCheck()
                functionspace=self.Ensemble[0][0][0].function_space()
                zerofunction=Function(functionspace)
                NewEnsemble=[]
                for i in range(len(self.Ensemble)):
                    N=len(self.Ensemble[i]) # get sample size
                    R=list(np.shape(zerofunction.dat.data)) # find dimension of data 
                    R.append(N)
                    GridDimensions=tuple(R)
                    GridC=np.zeros(GridDimensions)
                    GridF=np.zeros(GridDimensions)
                    for j in range(len(self.Ensemble[i])):
                        GridC[...,j]=self.Ensemble[i][j][0].dat.data
                        GridF[...,j]=self.Ensemble[i][j][1].dat.data
                    NewEnsemble.append(tuple([GridC,GridF]))
                self.Ensemble=NewEnsemble
                # turn on function indicator
                self.Type="Data"
            if To=="Function":
                NewEnsemble=[]
                for i in range(len(self.Ensemble)):
                    NewEnsemble.append([])
                    for j in range(np.shape(self.Ensemble[i][0])[-1]):
                        Fc=Function(self.__OriginalFunctionSpaces[i][0])
                        Ff=Function(self.__OriginalFunctionSpaces[i][1])
                        Fc.dat.data[:]=self.Ensemble[i][0][...,j]
                        Ff.dat.data[:]=self.Ensemble[i][1][...,j]
                        NewEnsemble[i].append(tuple([Fc,Ff]))
                self.Ensemble=NewEnsemble
                self.Type="Function"
        
        
class CopyEnsembleHierarchy():
    
    def __init__(self,Ensemble_Hierarchy):
        self.Copy=EnsembleHierarchy()
        self.Copy.__dict__.update(Ensemble_Hierarchy.__dict__)









