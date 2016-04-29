"""


State Variable Tuples


Alastair Gregory

"""


from packages import *


class state():
    def __init__(self,input_1,input_2):
        self.state=tuple([input_1,input_2])
        # give the state the attributes the levels of each fine / coarse solution. be careful for lists of states
        if type(self.state) is list:
            self.levels=tuple([get_level(self.state[0][0])[1],get_level(self.state[1][0])[1]])
        else:
            self.levels=tuple([get_level(self.state[0])[1],get_level(self.state[1])[1]])



