
from packages import *

''' Need to make ErrorBounding() take Firedrake Functions '''



def ErrorBoundingAttributeChecks(EnsembleHierarchy):
    E=EnsembleHierarchy
    from firedrake import Function
    if hasattr(E,'Mean')==0:
        raise AttributeError('EnsembleHierarchy Does Not Have Sample Statistics As Attribute')
    if hasattr(E,'Variance')==0:
        raise AttributeError('EnsembleHierarchy Does Not Have Sample Statistics As Attribute')
    if E.Type=='Data':
        if type(E.Mean[0])==Function:
            raise ValueError('Update SampleStatistics for new EnsembleHierarchy.Type')
    if E.Type=='Function':
        if type(E.Mean[0])!=Function:
            raise ValueError('Update SampleStatistics for new EnsembleHierarchy.Type')



def OptimalNl(EnsembleHierarchy,Eps,gamma=1):
    E=EnsembleHierarchy
    ErrorBoundingAttributeChecks(E)
    if E.Type=="Data":
        Tl=[]
        for i in range(E.L):
            Tl.append(np.linalg.norm(E.Variance[i]))
    if E.Type=="Function":
        Tl=[]
        for i in range(E.L):
            Tl.append(norm(E.Variance[i]))
    HL=np.array(E.hl)
    VL=np.array(Tl)
    Nl=np.ceil(np.multiply(np.sqrt(np.multiply(VL,HL**gamma)),((2*np.sum(np.sqrt(np.divide(VL,HL**gamma))))/(Eps**2)))) 
    return Nl



def Convergence(EnsembleHierarchy,Eps):
    E=EnsembleHierarchy
    ErrorBoundingAttributeChecks(E)
    if E.Type=="Data":
        EstBias=np.linalg.norm(E.Mean[-1])
    if E.Type=="Function":
        EstBias=norm(E.Mean[-1])
    Con=EstBias<((E.M-1)*Eps)/np.sqrt(2)
    return Con



class SampleStatistics():
    """ This is a wrapper to update sample state of ensemble hierarchy type """
    def __init__(self,EnsembleHierarchy): # This is wrapper for EnsembleHierarchy type structure
        self.__E = EnsembleHierarchy # inherit
        self.__Ensemble=self.__E.Ensemble
        # need to check which type the ensemble hierarchy is "Data" or "Function"
        ''' Update sample stats '''
        self.__E.Mean=self.__OnlineMean()
        self.__E.Variance=self.__OnlineCovariance()
        self.__E.MultilevelExpectation=self.__MultilevelMean()
    def __MultilevelMean(self):
        L=len(self.__Ensemble)
        if self.__E.Type=="Function":
            MLMCMean=Function(self.__Ensemble[0][0][0].function_space()).assign(0)
            for i in range(L):
                if i==0:
                    for j in range(len(self.__Ensemble[0])):
                        MLMCMean+=(1.0/float(len(self.__Ensemble[0])))*self.__Ensemble[0][j][1]
                else:
                    for j in range(len(self.__Ensemble[i])):    
                        MLMCMean+=(1.0/float(len(self.__Ensemble[i])))*(self.__Ensemble[i][j][1]-self.__Ensemble[i][j][0])
            MM=assemble(MLMCMean)
        if self.__E.Type=="Data":
            MLMCMean=np.zeros(np.shape(self.__Ensemble[0][1][...,-1]))
            for i in range(L):
                if i==0:
                    MLMCMean+=np.mean(self.__Ensemble[0][1],axis=-1)
                else:
                    MLMCMean+=np.mean(self.__Ensemble[i][1]-self.__Ensemble[i][0],axis=-1)
            MM=MLMCMean
        return MM
    def __OnlineMean(self):
        if self.__E.Type=="Data":
            El=[]
            for i in range(len(self.__Ensemble)):
                if i==0:
                    El.append(np.mean(self.__Ensemble[i][1],axis=-1))
                else:
                    El.append(np.mean(self.__Ensemble[i][1]-self.__Ensemble[i][0],axis=-1))
        if self.__E.Type=="Function":
            El=[]
            for i in range(len(self.__Ensemble)):
                if i==0:
                    G=Function(self.__Ensemble[i][0][1].function_space()).assign(0)
                    for j in range(len(self.__Ensemble[i])):
                        G+=(1/float(len(self.__Ensemble[i])))*self.__Ensemble[i][j][1]
                    El.append(assemble(G))
                else:
                    Gc=Function(self.__Ensemble[i][0][0].function_space()).assign(0)
                    Gf=Function(self.__Ensemble[i][0][1].function_space()).assign(0)
                    for j in range(len(self.__Ensemble[i])):
                        Gf+=(1/float(len(self.__Ensemble[i])))*(self.__Ensemble[i][j][1])
                        Gc+=(1/float(len(self.__Ensemble[i])))*(self.__Ensemble[i][j][0])
                    El.append(assemble(Gf-Gc))
        return El
    """ Used For Consistency Checks """ """
    def OnlineGoldStandardMean(self,FunctionEnsemble):
        # this one is really only for confirming consistency, can get rid of!
        m=0
        for i in range(len(FunctionEnsemble[-1])):
            m+=FunctionEnsemble[-1][i][1].at([0.5,0.5])*(1/len(FunctionEnsemble[-1]))
        return m
    def OnlineEnsembleMean(self,FunctionEnsemble):
        # this one is really only for confirming consistency, can get rid of!
        m=np.zeros(len(FunctionEnsemble))
        for i in range(len(FunctionEnsemble)):
            for j in range(len(FunctionEnsemble[i])):
                if i==0:
                    m[i]+=FunctionEnsemble[i][j][1].at([0.5,0.5])*(1/len(FunctionEnsemble[i]))
                else:
                    m[i]+=(FunctionEnsemble[i][j][1].at([0.5,0.5])-FunctionEnsemble[i][j][0].at([0.5,0.5]))*(1/len(FunctionEnsemble[i]))
        return m
    """ """    *************************   """
    def __OnlineCovariance(self): # THIS NEEDS TO BE CHANGED TO ACTUALLY BE COVARIANCE -> NOT JUST DIAGONAL OF VAR OF ALL COMPONENTS -> DO SOMETHING LIKE IN SPATIAL POST PROCESSING
        if self.__E.Type=="Data":
            Vl=[]
            for i in range(len(self.__Ensemble)):
                N=np.shape(self.__Ensemble[i][0])[-1]
                e0=np.ravel(np.mean(np.square(self.__Ensemble[i][1]-self.__Ensemble[i][0]),axis=-1)-np.square(np.mean(self.__Ensemble[i][1]-self.__Ensemble[i][0],axis=-1)))
                # e0 has been changed for actual difference between variances!!
                e1=np.ravel(np.mean(np.square(self.__Ensemble[i][1]),axis=-1)-np.square(np.mean(self.__Ensemble[i][1],axis=-1)))
                if i==0:
                    Vl.append(e1)
                else:
                    Vl.append(e0)
        if self.__E.Type=="Function":
            Vl=[]
            for i in range(len(self.__Ensemble)):
                if i==0:
                    G=Function(self.__Ensemble[i][0][1].function_space()).assign(0)
                    GSq=Function(self.__Ensemble[i][0][1].function_space()).assign(0)
                    for j in range(len(self.__Ensemble[i])):
                        G+=(1/float(len(self.__Ensemble[i])))*self.__Ensemble[i][j][1]
                        GSq+=(1/float(len(self.__Ensemble[i])))*(self.__Ensemble[i][j][1]**2)
                    Vl.append(assemble(GSq-(G**2)))  # append actual functions
                else:
                    Gc=Function(self.__Ensemble[i][0][0].function_space()).assign(0)
                    Gf=Function(self.__Ensemble[i][0][1].function_space()).assign(0)
                    GSq=Function(self.__Ensemble[i][0][1].function_space()).assign(0)
                    for j in range(len(self.__Ensemble[i])):
                        Gf+=(1/float(len(self.__Ensemble[i])))*(self.__Ensemble[i][j][1])
                        Gc+=(1/float(len(self.__Ensemble[i])))*(self.__Ensemble[i][j][0])
                        GSq+=(1/float(len(self.__Ensemble[i])))*((self.__Ensemble[i][j][1]-self.__Ensemble[i][j][0])**2)
                    Vl.append(assemble(GSq-((Gf-Gc)**2))) # append actual functions. Here sample variance \in V_{l} <- the function space that the level to append to is on
        return Vl
    





