"""

Prolong and Injecting

Alastair Gregory 2016

"""
from packages import *


class ProlongInject():
    """ Allows one to access functions to prolong / inject over mutiple levels, with hierarchal function inputs """
    def TurnHierarchyIntoItemList(self,Hierarchy):
        """
        """
        L=len(Hierarchy)
        List=[]
        for i in range(L):
            List.append(Function(Hierarchy[i].function_space()))
        return List
    def ProlongUpToFinestLevel(self,CoarseIndex,Function,Hierarchy):
        """
        """
        Hier=self.TurnHierarchyIntoItemList(Hierarchy)
        Hier[CoarseIndex]=Function
        if CoarseIndex==len(Hierarchy):
            Hier[-1]=Function
        for i in range(len(Hier)-1-CoarseIndex):
            prolong(Hier[i+CoarseIndex],Hier[i+1+CoarseIndex])
        return Hier[-1]
    def ProlongUpToAnyLevel(self,CoarseIndex,Level_to_prolong_to,Function,Hierarchy): 
        """
        """
        if CoarseIndex>Level_to_prolong_to:
            raise ValueError('Cant prolong to level less than current')
        Hier=self.TurnHierarchyIntoItemList(Hierarchy)
        Hier[CoarseIndex]=Function
        if CoarseIndex==Level_to_prolong_to: # this is the same as Level_to_prolong_to
            Hier[Level_to_prolong_to]=Function
        for i in range(Level_to_prolong_to-CoarseIndex):
            prolong(Hier[CoarseIndex+i],Hier[CoarseIndex+i+1])
        return Hier[Level_to_prolong_to]
    def InjectDownToAnyLevel(self,FineIndex,Level_to_inject_to,Function,Hierarchy): 
        """
        """
        if FineIndex<Level_to_inject_to:
            raise ValueError('Cant inject to level greater than current')
        Hier=self.TurnHierarchyIntoItemList(Hierarchy)
        Hier[FineIndex]=Function
        if FineIndex==Level_to_inject_to: 
            Hier[Level_to_inject_to]=Function
        for i in range(FineIndex-Level_to_inject_to):
            inject(Hier[FineIndex-i],Hier[FineIndex-i-1])
        return Hier[Level_to_inject_to]


    
    
        
        
        




