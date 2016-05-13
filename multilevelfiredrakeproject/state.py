
from packages import *


class state():
    
    """Tuple that stores the state :class:`Function` and properties of them for the coarse (index 0) and fine (index 1) discretizations for a single realisation of the user define system.
    
    	:param input_1: Coarse state :class:`Function` needed for discretization. Can be a :attr:`list` of multiple :class:`Function`.
    	:param input_2: Fine state :class:`Function` needed for discretization. Can be a :attr:`list` of multiple :class:`Function`.
    	:type input_1: :class:`Function`.
    	:type input_2: :class:`Function`.
    
    """
    
    def __init__(self,input_1,input_2):
        self.state=tuple([input_1,input_2]) #: :attr:`tuple` of coarse and fine state :class:`Function`
        # give the state the attributes the levels of each fine / coarse solution. be careful for lists of states.
        # preallocate prepared state
        self.prepared_state=None
        if type(self.state) is list:
            self.levels=tuple([get_level(self.state[0][0])[1],get_level(self.state[1][0])[1]]) #: the levels corresponding to the coarse and fine state :class:`Function`.
        else:
            self.levels=tuple([get_level(self.state[0])[1],get_level(self.state[1])[1]])
        # add check for levels - non fatal!!
        if self.levels[0]==-1 or self.levels[1]==-1:
            raise Warning('Levels of state may not be actual hierarchal levels. Check if they belong to FunctionHierarchy! get_level has failed.')
        



