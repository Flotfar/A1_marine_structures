

from Input_units_A1 import *
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

#### Migration speeds ####
d = d_50
V_h1, V_h2 = migration_speed_span_shoulders(regime, U_f, nu, g, s, d, e, D)
print("Primary migration speed = " + str(round(V_h1, 3)) + "\n" + "Secondary migration speed = " + str(round(V_h2, 3)) + "\n")

#Time from primary to secondary

