

from packages import *



class EnsembleHierarchy():
    def __init__(self,M=2):
        self.Type="Function"
        self.Ensemble=[]
        self.hl=[]
        self.L=0
        self.M=M
        # Preallocate attributes
        self.MultilevelExpectation=None
        self.Mean=[]
        self.Variance=[]
        #
        self.__OriginalFunctionSpaces=[]
        self.Weights=[] # added weights into this (evenly weighted is default! Will be easy to add importance sampling when it's added in!
    def __SameFunctionSpaceErrorCheck(self):
        """" Checks for the same function space for every ensemble member in hierarchy """
        Comparison=len(self.Ensemble[0][0][0].dat.data)
        for i in range(len(self.Ensemble)):
            for j in range(len(self.Ensemble[i])):
                for k in range(len(self.Ensemble[i][j])):
                    if len(self.Ensemble[i][j][k].dat.data)!=Comparison:
                        raise ValueError('All Ensemble Members Do Not Have Same Function Spaces')
    def AppendToEnsemble(self,tuple_to_append,Level_to_append_to):
        """ Append into the Level_to_append_to index, if exists, of ensemble hierarchy
        Inputs:
        - Ensemble, LIST
        - tuple_to_append, tuple to append to Ensemble in level...
        - Level_to_append_to, level, if exists, to append to Ensemble
        Outputs:
        - UpdatedEnsemble
        """
        if hasattr(tuple_to_append,'prepared_state')==0:
            raise AttributeError('State hasnt been prepared for appedning into ensemble hierarchy')
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
        """ Transfers Ensemble Into Structured Grid of Cell Data """
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









