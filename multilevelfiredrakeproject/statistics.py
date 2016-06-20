
from packages import *

''' Need to make ErrorBounding() take Firedrake Functions '''


def ErrorBoundingAttributeChecks(EnsembleHierarchy):
    
    """

    """

    E = EnsembleHierarchy
    
    from firedrake import Function
    
    if hasattr(E, 'Mean') == 0:
        raise AttributeError(
            'EnsembleHierarchy Does Not Have Sample Statistics As Attribute')
    
    if hasattr(E, 'Variance') == 0:
        raise AttributeError(
            'EnsembleHierarchy Does Not Have Sample Statistics As Attribute')
    
    if not isinstance(E.Mean[0], Function):
        raise ValueError(
            'Update SampleStatistics for new EnsembleHierarchy.Type')


def OptimalNl(EnsembleHierarchy, Eps, gamma=1):
    
    """

    """

    E = EnsembleHierarchy
    
    ErrorBoundingAttributeChecks(E)
    
    Tl = []
    
    for i in range(E.L):
        Tl.append(norm(E.Variance[i]))
    HL = np.array(E.nxl)**(-1)
    VL = np.array(Tl)
    
    Nl = np.ceil(np.multiply(np.sqrt(np.multiply(VL, HL**gamma)),
                             ((2 * np.sum(np.sqrt(np.divide(VL, HL**gamma)))) / (Eps**2))))
    return Nl


def Convergence(EnsembleHierarchy, Eps):
    
    """

    """

    
    E = EnsembleHierarchy
    
    ErrorBoundingAttributeChecks(E)
    
    EstBias = norm(E.Mean[-1])
    Con = EstBias < ((E.M - 1) * Eps) / np.sqrt(2)
    
    return Con


class SampleStatistics(object):

    """ This is a wrapper to update sample state of ensemble hierarchy type

    """

    # This is wrapper for EnsembleHierarchy type structure
    def __init__(self, EnsembleHierarchy):
        
        """

        """

        self.__E = EnsembleHierarchy  # inherit
        
        ''' Update sample stats every time class is initialized '''
        
        c = 0
        if self.__E.Type == 'Data':
            self.__E.EnsembleTransfer('Function')
            # transfer to get sample stats within form of functions
            c = 1
        
        self.__Ensemble = self.__E.Ensemble
        
        self.__E.Mean = self.__OnlineMean()
        self.__E.Variance = self.__OnlineCovariance()
        self.__E.MultilevelExpectation = self.__MultilevelMean()
        
        if c == 1:
            self.__E.EnsembleTransfer('Data')  # transfer back

        super(SampleStatistics, self).__init__()

    def __MultilevelMean(self):
        
        """

        """

        L = len(self.__Ensemble)
        
        if self.__E.Type == "Function":
            MLMCMean = Function(
                self.__Ensemble[0][0][0].function_space()).assign(0)
            for i in range(L):
                
                if i == 0:
                    for j in range(len(self.__Ensemble[0])):
                        MLMCMean += (1.0 /
                                     float(len(self.__Ensemble[0]))) * self.__Ensemble[0][j][1]
                
                else:
                    for j in range(len(self.__Ensemble[i])):
                        MLMCMean += (1.0 / float(len(self.__Ensemble[i]))) * (
                            self.__Ensemble[i][j][1] - self.__Ensemble[i][j][0])
            
            MM = assemble(MLMCMean)
        
        else:
            
            raise TypeError(
                'ensemble hierarchy has not been temporarily transfered to function type')
        
        return MM

    def __OnlineMean(self):
        
        """

        """

        if self.__E.Type == "Function":
            
            El = []
            for i in range(len(self.__Ensemble)):
                
                if i == 0:
                    
                    G = Function(
                        self.__Ensemble[i][0][1].function_space()).assign(0)
                    
                    for j in range(len(self.__Ensemble[i])):
                        G += (1 /
                              float(len(self.__Ensemble[i]))) * self.__Ensemble[i][j][1]
                    
                    El.append(assemble(G))
                
                else:
                    
                    Gc = Function(
                        self.__Ensemble[i][0][0].function_space()).assign(0)
                    Gf = Function(
                        self.__Ensemble[i][0][1].function_space()).assign(0)
                    
                    for j in range(len(self.__Ensemble[i])):
                        Gf += (1 / float(len(self.__Ensemble[i]))
                               ) * (self.__Ensemble[i][j][1])
                        Gc += (1 / float(len(self.__Ensemble[i]))
                               ) * (self.__Ensemble[i][j][0])
                    
                    El.append(assemble(Gf - Gc))
        
        else:
            
            raise TypeError(
                'ensemble hierarchy has not been temporarily transfered to function type')
        
        return El

    # THIS NEEDS TO BE CHANGED TO ACTUALLY BE COVARIANCE -> NOT JUST DIAGONAL
    # OF VAR OF ALL COMPONENTS -> DO SOMETHING LIKE IN SPATIAL POST PROCESSING
    def __OnlineCovariance(self):
        
        """


        """

        if self.__E.Type == "Function":
            
            Vl = []
            for i in range(len(self.__Ensemble)):
                if i == 0:
                    
                    G = Function(
                        self.__Ensemble[i][0][1].function_space()).assign(0)
                    
                    GSq = Function(
                        self.__Ensemble[i][0][1].function_space()).assign(0)
                    
                    for j in range(len(self.__Ensemble[i])):
                        G += (1 /
                              float(len(self.__Ensemble[i]))) * self.__Ensemble[i][j][1]
                        
                        GSq += (1 / float(len(self.__Ensemble[i]))) * (
                            self.__Ensemble[i][j][1]**2)
                    
                    # append actual functions
                    Vl.append(assemble(GSq - (G**2)))
                
                else:
                    
                    Gc = Function(
                        self.__Ensemble[i][0][0].function_space()).assign(0)
                    Gf = Function(
                        self.__Ensemble[i][0][1].function_space()).assign(0)
                    
                    GSq = Function(
                        self.__Ensemble[i][0][1].function_space()).assign(0)
                    
                    for j in range(len(self.__Ensemble[i])):
                        Gf += (1 / float(len(self.__Ensemble[i]))
                               ) * (self.__Ensemble[i][j][1])
                        Gc += (1 / float(len(self.__Ensemble[i]))
                               ) * (self.__Ensemble[i][j][0])
                        
                        GSq += (1 / float(len(self.__Ensemble[i]))) * (
                            (self.__Ensemble[i][j][1] - self.__Ensemble[i][j][0])**2)
                    
                    # append actual functions. Here sample variance \in V_{l}
                    # <- the function space that the level to append to is on
                    Vl.append(assemble(GSq - ((Gf - Gc)**2)))
        
        else:
            
            raise TypeError(
                'ensemble hierarchy has not been temporarily transfered to function type')
        
        return Vl
