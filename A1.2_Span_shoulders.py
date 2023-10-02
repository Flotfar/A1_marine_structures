

from A1_pipe_Input_units import *
from General_constants import *
from Scour_functions import *

####################################################
#### Input Units
####################################################

#### General units ####
kappa, ks, g, nu, U_f, s, H_s, Fr = general_const(d_50, V, rho_s, rho_w, h, H_s, H)

#### Case specific units ####
regime = "live_bed"    #either "live_bed" or "clear_water"


####################################################
#### Results
####################################################
print("\n######## Results: ########")

#### Migration speeds and time from primary to secondary ####
d = d_50
V_h1, V_h2, time = migration_span_shoulders(regime, U_f, nu, g, s, d, e, D)
