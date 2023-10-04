import numpy as np

from A6_Input_units import *
from A1_6_Monopile_design import *
from General_constants import *
from Scour_functions import *
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#From fig 12 in Nielsen and Petersen, find the critical shield parameters:

#summer
theta_crs = 1.5*10**(-2)

#design:
theta_crd = 3.1*10**(-2)

#Critical theta will then be:
theta_cr = min(theta_crs,theta_crd)

#Summer is critical.

#b)

##############
## Variables##
##############

Wb = 4*D
hb = 1 #m 

#Stone width
Dc = 0.103 #m

# Stone layers 
N = hb/Dc

#Berm height to hberm width aspect ratio (4) in Petersen et al 2015
Ar = hb/Wb

#Friction velocity from design conditions
U_f = Ud/(6+2.5*np.log(h/(2.5*d50)))

#Shield parameter based on friction velocity from design conditions is:

theta = U_f**2/(g*(s-1)*d50)
