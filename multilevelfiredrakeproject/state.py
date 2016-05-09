"""


State Variable Tuples


Alastair Gregory

"""


from packages import *


class state():
    
    """Tuple that stores the state :class:`Functions' and properties of them for the coarse (index 0) and fine (index 1) discretizations for a single realisation of the user define system
    
    	:arg input_1: Coarse state :class:`Functions' needed for discretization. Can be a list of :class:`Functions'.
    	:arg input_2: Fine state :class:`Functions' needed for discretization. Can be a list of :class:`Functions'.
    
    """
    
    def __init__(self,input_1,input_2):
        self.state=tuple([input_1,input_2])
        # give the state the attributes the levels of each fine / coarse solution. be careful for lists of states
        if type(self.state) is list:
            self.levels=tuple([get_level(self.state[0][0])[1],get_level(self.state[1][0])[1]])
        else:
            self.levels=tuple([get_level(self.state[0])[1],get_level(self.state[1])[1]])
        # add check for levels - non fatal!!
        if self.levels[0]==-1 or self.levels[1]==-1:
            raise Warning('Levels of state may not be actual hierarchal levels. Check if they belong to FunctionHierarchy! get_level has failed.')
        



