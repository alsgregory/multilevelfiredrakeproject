"""


Ensemble Forecasting


Alastair Gregory

"""
from packages import *

class EnsembleForecast():
    def __init__(self,EnsembleHierarchy,Weights,DesiredFunctionSpace,N,localisation=False,Sigma=1):
        self.__OriginalFunctionSpace = EnsembleHierarchy._EnsembleHierarchy__OriginalFunctionSpaces[-1][1]
        # Needs to be in Data Type
        if EnsembleHierarchy.Type!="Data":
            EnsembleHierarchy.EnsembleTransfer(To="Data")
            transfered=1
        else:
            transfered=0
        self.requiredshape=np.shape(EnsembleHierarchy.Ensemble[-1][1][...,-1])
        #
        NE=[]
        for i in range(len(EnsembleHierarchy.Ensemble)):
            NE.append([])
            for j in range(len(EnsembleHierarchy.Ensemble[i])):
                NE[i].append(self.__FlattenMultiD(EnsembleHierarchy.Ensemble[i][j]))
        self.Ensemble=NE
        # Define forecast ensemble size
        self.N=N
        self.Weights=Weights
        self.dimension_of_state_space=self.__OriginalFunctionSpace.mesh().geometric_dimension()
        self.localisation=localisation
        ''' Need OriginalFunctionSpace to be DG if want to use localisation OR NON-VECTOR!! '''
        if self.localisation==True:
            if self.__OriginalFunctionSpace.ufl_element().family()=='Lagrange' or self.__OriginalFunctionSpace.ufl_element().value_size()>1:
                raise ValueError('Cannot use localisation with CG or Vector elements. Please convert to DG (and non-vector) and change DesiredFunctionSpace to CG (or vector) as to revert back')
        if self.localisation==True:
            self.Grid=self.__GenerateDistanceGrid()
        else:
            self.Grid=None
        self.Sigma=Sigma
        # Preallocate Delta attribute
        self.Delta=None
        self.Cov=self.__MLMCCovariance(self.Ensemble) # define covariance matrix for object - This may not be pos def, so use shrinking!
        #self.c=np.cov(self.Ensembles[-1][1])
        self.__S=SpatialPostProcessing(self.Grid,self.dimension_of_state_space,self.Cov,self.N,self.localisation,self.Sigma)
        self.DesiredFunctionSpace = DesiredFunctionSpace
        self.Forecast=self.__GenerateEnsemble()
        self.Type='Function' # Default type after generated enssemble
        # convert back to function
        if transfered==1:
            EnsembleHierarchy.EnsembleTransfer(To="Function")
        # revert ensemble
        AE=[]
        for i in range(len(self.Ensemble)):
            AE.append([])
            for j in range(len(self.Ensemble[i])):
                AE[i].append(self.__ReturnMultiD(self.Ensemble[i][j],self.requiredshape))
        self.Ensemble=AE
    def __WeightedSampleQuantile(self,u,X,W):
        """ Finds the consistent approximation to the quantile function F^(-1)(u)=x where u~U[0,1], x~X, for given u 
        Inputs:
        - u, uniform sample, numpy array
        - X, numpy array sized N
        Ouputs:
        - x, numpy array of approximations to F^(-1)(u)
        """
        A=np.argsort(X)
        x=np.sort(X)
        w=W[A.astype(int)]
        CW=np.cumsum(w)
        i=np.searchsorted(CW,u)
        return x[i.astype(int)]
    def __IterativeDistance(self,U,f,index):
		# project into DG0
		old_degree=f.function_space().ufl_element().degree()
		F=Function(FunctionSpace(U,'DG',0)).project(f)
		YCG=FunctionSpace(U,'CG',2)
		#
		cellaverage2vertex_kernel = """
		for(int i=0;i<vertex.dofs;i++){
		vertex[i][0]+=cell[0][0];
		}
		"""
		vertex2cellaverage_kernel = """ float scale=0.0;
		for(int i=0;i<vertex.dofs;i++){
		new_cell[0][0]+=vertex[i][0];
		scale=scale+1.0;
		}
		new_cell[0][0]=new_cell[0][0]/scale;
		"""
		#
		n=len(F.dat.data)
		ADG=FunctionSpace(U,'DG',0)
		i=0
		while np.any(F.dat.data==0.0):
		    vertex_cell=Function(YCG)
		    vertex_cell.assign(0.0)
		    par_loop(cellaverage2vertex_kernel,dx,{"vertex":(vertex_cell,RW),"cell":(F,READ)})
		    new_cell=Function(ADG)
		    new_cell.assign(0.0)
		    par_loop(vertex2cellaverage_kernel,dx,{"new_cell":(new_cell,RW),"vertex":(vertex_cell,READ)})
		    F.assign(new_cell)
		    i+=1
		F_new=Function(FunctionSpace(U,'DG',old_degree)).project(F)
		distance=np.ones(len(F_new.dat.data))-np.divide(F_new.dat.data,F_new.dat.data[index]*np.ones(len(F_new.dat.data)))
		distance=np.multiply(distance,(distance>=0))
		return distance
    def __GenerateDistanceGrid(self):
        ''' Only in DG and non-vector '''
        # mesh
        U=self.__OriginalFunctionSpace.mesh()
        n=len(self.Ensemble[-1][1][:,0])
        distance=np.zeros((n,n))
        for i in range(n):
            f=Function(self.__OriginalFunctionSpace)
            f.dat.data[i]=1.0
            distance[i,:]=self.__IterativeDistance(U,f,i)
        return distance
    def EnsembleTransfer(self,To):
        if To==self.Type:
            ''' nothing happens '''
        if To!=self.Type:
            if To=='Data':
                N=len(self.Forecast)
                NewEnsemble=np.zeros((len(self.Forecast[0].dat.data),N))
                for i in range(N):
                    NewEnsemble[:,i]=self.Forecast[i].dat.data
                self.Forecast=NewEnsemble
                self.Type='Data'
            if To=='Function':
                self.Forecast=self.__RevertEnsemble(self.Forecast)
                self.Type='Function'
    def __MLMCCovariance(self,Ensembles):
        for i in range(len(Ensembles)):
            if i==0:
                C=np.cov(Ensembles[0][1])
            else:
                C+=np.cov(Ensembles[i][1])-np.cov(Ensembles[i][0])
        return C
    def __FlattenMultiD(self,MultiDArray):
        """
        Turns a greater than 2 dimension (D x ... x D x N) numpy array into a (D x N) flattened array
        """
        FlattenedArray=np.zeros((np.prod(np.shape(MultiDArray[...,0])),np.shape(MultiDArray)[-1]))
        for i in range(np.shape(MultiDArray)[-1]):
            FlattenedArray[:,i]=np.ravel(MultiDArray[...,i],order='C')
        return FlattenedArray
    def __RevertEnsemble(self,GridEnsemble):
        """ Transfers an ensemble forecast D x N back into functions, given the level of the hierarchy it's on """
        FunctionForecasts=[] # LIST
        for i in range(len(GridEnsemble[0,:])):
            # revert back to old fs.
            F=Function(self.__OriginalFunctionSpace)
            F.dat.data[...]=GridEnsemble[...,i]
            # now project to desired function space
            FNew=Function(self.DesiredFunctionSpace)
            FNew.project(F)
            FunctionForecasts.append(FNew)
        return FunctionForecasts
    def __ReturnMultiD(self,FlattenedArray,requiredshape):
        """
        Returns a flattened array back to multi d
        Inputs:
        - FlattenedArray (D x N)
        - requiredshape tuple>=2
        """
        l=list(requiredshape)
        l.append(len(FlattenedArray[0,:]))
        R=tuple(l)
        Array=np.zeros(R)
        for i in range(np.shape(FlattenedArray)[-1]):
            Array[...,i]=np.reshape(FlattenedArray[:,i],(requiredshape),order='C')
        return Array
    def __GenerateEnsemble(self): 
        n=len(self.Ensemble[0][0][0,:]) # size of new single ensemble (desired)
        L=len(self.Ensemble)
        D=len(self.Ensemble[0][0][:,0])
        Uniforms=np.random.uniform(0,1,((D,self.N))) # produce random samples using ITS then correlate
        print colored('Generating a multilevel ensemble forecast, length '+str(self.N)+'...','red')
        for d in range(D):
            Uniforms[d,:]=np.linspace(0,self.N-1,self.N)/self.N
            np.random.shuffle(Uniforms[d,:])
        # First sample uncorrelated marginals
        UncorrelatedSamples=np.zeros((D,self.N))
        for d in range(D):
            for i in range(L):
                if i==0:
                    UncorrelatedSamples[d,:]+=self.__WeightedSampleQuantile(Uniforms[d,:],self.Ensemble[0][1][d,:],self.Weights[0][1][:])
                else:
                    UncorrelatedSamples[d,:]+=(self.__WeightedSampleQuantile(Uniforms[d,:],self.Ensemble[i][1][d,:],self.Weights[i][1][:])-self.__WeightedSampleQuantile(Uniforms[d,:],self.Ensemble[i][0][d,:],self.Weights[i][0][:]))
        # Generate Copula Uniform Samples
        CopulaSamples=self.__S._SpatialPostProcessing__JointNormalTransform()
        # Update covariance matrix
        self.Cov=self.__S.C
        # Update delta used
        self.Delta=self.__S.Delta
        # if dimension is greater than 1, copula the uncorrelated sample
        CorrelatedSamples=np.zeros((D,self.N))
        for d in range(D):
            CorrelatedSamples[d,:]=self.__WeightedSampleQuantile(CopulaSamples[d,:],UncorrelatedSamples[d,:],np.ones(self.N)*(1.0/float(self.N)))
        print colored('****** Ensemble Forecasted ******','red')
        return self.__RevertEnsemble(self.__ReturnMultiD(CorrelatedSamples,self.requiredshape))



class SpatialPostProcessing():
    def __init__(self,Grid,Dimension,CovarianceMatrix,N,T=False,Sigma=1):
        self.dimension_of_state_space=Dimension
        self.Grid=Grid
        ''' This grid is either a tuple of 2d grids or an array 1 grid '''
        self.C=CovarianceMatrix
        self.Sigma=Sigma
        self.T=T
        self.N=N
    def __Tapering(self,C):
        """ Tapering method - Bickel et al 2008 """
        # Define distance matrix, only works for 1 and 2 dimension of state space.
        Distance=self.Grid
        # Checks if one needs to taper
        if self.T==True:
            R=np.exp(-Distance/self.Sigma)
            Regularized=C*R
        elif self.T==False:
            Regularized=C
        else:
            raise ValueError('Option to taper needs to be Boolean type')
        return Regularized
    def __Shrink(self,C,delta):
        """ Shrinking of an approximation to a correlation matrix, to guarantee that it is positive definite (Nick Highman 2014). Target is the identity
        Inputs:
        - C, Correlation Matrix
        - delta, Found by the bisection method
        Ouputs:
        - Shrunk correlation matrix, positive definite
        """
        return (delta*C)+((1-delta)*np.identity(len(C[0,:])))
    def __bisectionmethod(self,C):
        """ Finds optimal delta, that allows C to be positive definite under shrinking (Nick Highman 2014)
        Inputs:
        - C, Correlation Matrix
        Outputs:
        - Delta, to be used in Shrinking Function
        """
        al=0.0; ar=1.0; epsilon=1e-6
        if np.all(np.linalg.eigh(self.__Shrink(C,1))[0]>0)==1: # initial check
            return 1
        while ar-al>epsilon: # mid point method
            am=(ar+al)/2.0
            if np.all(np.linalg.eigh(self.__Shrink(C,am))[0]>0)==1:
                al=am
            else:
                ar=am
        return al
    def __CovintoCorr(self,Cov):
        D=len(Cov[0,:])
        Corr=np.zeros((D,D))
        for i in range(D):
            Corr[i,:]=np.divide(Cov[i,:],np.sqrt(np.multiply(Cov[i,i],np.diag(Cov))))
        return Corr
    def __JointNormalTransform(self):
        """ This samples N Uniform points to invert marginals! """
        # convert into standard normal - correlation matrix!
        print colored('Finding delta needed to shrink mutlilevel covariance...')
        # Carry out localisation and shrinking to guarantee Pos Def - Print Delta needed!
        delta=self.__bisectionmethod(self.C)
        print colored('Delta Needed To Shrink MLMC Cov: '+str(delta),'green')
        self.Delta=delta # save delta used for forecast
        print colored('Shrinking multilevel covariance...','green')
        self.C=self.__Shrink(self.C,delta)
        if self.T==True:
            print colored('Tapering multilevel covariance...','green')
            self.C=self.__Tapering(self.C)
            print colored('Tapered and Shrunk','green')
        else:
            print colored('Shrunk','green')
        # Sample from MVN
        Samples=np.random.multivariate_normal(np.zeros(len(self.C[0,:])),self.C,self.N).T
        # CDF them
        print colored('Finiding Gaussian copula samples...','green')
        Uniforms=np.zeros(np.shape(Samples))
        Sds=np.diag(self.C)
        for i in range(len(Samples[:,0])):
            Uniforms[i,:]=NORM.cdf(Samples[i,:],0,np.sqrt(Sds[i]))
        return Uniforms
    
    




