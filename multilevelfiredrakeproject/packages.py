from __future__ import division # Get proper divison
import numpy as np
import math
import random
from scipy import stats
from scipy.stats import norm as NORM
from firedrake import *
parameters["reorder_meshes"] = False
from firedrake.mg.utils import *
from firedrake.mg.utils import get_level
from colorama import init
init()
from colorama import Fore, Back, Style
from termcolor import colored, cprint
