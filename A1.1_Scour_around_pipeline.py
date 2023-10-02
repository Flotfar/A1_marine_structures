
import numpy as np

from A1_pipe_Input_units import *
from General_constants import *
from Scour_functions import *

####################################################
#### Input Units
####################################################

# General units
kappa, ks, g, nu, U_f, s, H_s, Fr = general_const(d_50, V, rho_s, rho_w, h, H_s, H)


####################################################
#### Results
####################################################

print("\n######## Results: ########")
#### Scour for a steady current ####
U, U_cr = pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g)

#### Scour from waves: ####
U_m, KC = pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e)

print("##########################\n")

# when gaining the result input the reading (input a pop up menu)
reading = 0.05
U_mcr = pipeline_scour_waves_Ucr(reading, g, D, n, s)
