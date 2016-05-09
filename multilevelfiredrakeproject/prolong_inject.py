"""

Prolong and Injecting

Alastair Gregory 2016

"""
from packages import *

# Make sure that we don't explicitely take as arguments the levels of these functions as we can just work these out from get_levels. Remember warnings about is they are part of hierarchy or not!!!! -> test this!


class ProlongInject():
    
    
    """ Contains methods that are used for iteratively prolonging and inject coarse / fine :class:`Functions' up more than one level"""
    
    
    def __TurnHierarchyIntoItemList(self,Hierarchy):
        
        L=len(Hierarchy)
        List=[]
        for i in range(L):
            List.append(Function(Hierarchy[i].function_space()))
        return List
    
    
    def ProlongUpToFinestLevel(self,CoarseIndex,Function,Hierarchy):
        
        """Prolongs any :class:`Function' in a :class:`FunctionHierarchy' from level :arg CoarseIndex: to the finest level in the empty output :class:`FunctionHierarchy'
        
        	:arg CoarseIndex:
        	:arg :class:`Function:
        	:arg :class:`FunctionHierarchy': To inset the prolonged :class:`Function' into
        
        """
        
        Hier=self.__TurnHierarchyIntoItemList(Hierarchy)
        Hier[CoarseIndex]=Function
        if CoarseIndex==len(Hierarchy):
            Hier[-1]=Function
        for i in range(len(Hier)-1-CoarseIndex):
            prolong(Hier[i+CoarseIndex],Hier[i+1+CoarseIndex])
        return Hier[-1]
    
    
    def ProlongUpToAnyLevel(self,CoarseIndex,Level_to_prolong_to,Function,Hierarchy): 
        
        """Prolongs any :class:`Function' in a :class:`FunctionHierarchy' from level :arg CoarseIndex: to the level :arg Level_to_prolong_to: in the empty output :class:`FunctionHierarchy'
        
        	:arg CoarseIndex:
        	:arg Level_to_prolong_to:
        	:arg :class:`Function:
        	:arg :class:`FunctionHierarchy': To inset the prolonged :class:`Function' into
        
        """
        
        if CoarseIndex>Level_to_prolong_to:
            raise ValueError('Cant prolong to level less than current')
        Hier=self.__TurnHierarchyIntoItemList(Hierarchy)
        Hier[CoarseIndex]=Function
        if CoarseIndex==Level_to_prolong_to: # this is the same as Level_to_prolong_to
            Hier[Level_to_prolong_to]=Function
        for i in range(Level_to_prolong_to-CoarseIndex):
            prolong(Hier[CoarseIndex+i],Hier[CoarseIndex+i+1])
        return Hier[Level_to_prolong_to]
    
    
    def InjectDownToAnyLevel(self,FineIndex,Level_to_inject_to,Function,Hierarchy): 
        
        """Injects any :class:`Function' in a :class:`FunctionHierarchy' from level :arg FineIndex: to the level :arg Level_to_inject_to: in the empty output :class:`FunctionHierarchy'
        
        	:arg FineIndex:
        	:arg Level_to_inject_to:
        	:arg :class:`Function:
        	:arg :class:`FunctionHierarchy': To inset the injected :class:`Function' into
        
        """
        
        if FineIndex<Level_to_inject_to:
            raise ValueError('Cant inject to level greater than current')
        Hier=self.__TurnHierarchyIntoItemList(Hierarchy)
        Hier[FineIndex]=Function
        if FineIndex==Level_to_inject_to: 
            Hier[Level_to_inject_to]=Function
        for i in range(FineIndex-Level_to_inject_to):
            inject(Hier[FineIndex-i],Hier[FineIndex-i-1])
        return Hier[Level_to_inject_to]


    
    
        
        
        




