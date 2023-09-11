
from A1_Input_units import *
from General_constants import *
from Scour_functions import *
from Plots import *

####################################################
#### Input Units
####################################################

#### General units ####
kappa, ks, g, nu, U_f, s, H_s, Fr = general_const(d_50, V, rho_s, rho_w, h, H_s, H)


####################################################
#### Results
####################################################

#### Current alone ####
condition = "Steady current"    # 1) "Steady current"  2) "Waves"  3) "Tidal current"  4) "Waves + Currents"
#Determineing condition parameters:
U, U_cr, U_c = pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g)
U_m = 0
KC = 0
#Determining the scour profile values:
S_eq, W1, W2 = scour_profile(condition, e, D, KC, U_c, U_m)

#plotting the scour profile
Scour_profile_plot(D, W1, W2, S_eq, condition)



#### Tidal current ####
condition = "Tidal current"    # 1) "Steady current"  2) "Waves"  3) "Tidal current"  4) "Waves + Currents"
#Determining the scour profile values:
S_eq, W1, W2 = scour_profile(condition, e, D, KC, U_c, U_m)
#plotting the scour profile
Scour_profile_plot(D, W1, W2, S_eq, condition)



#### Waves alone ####
condition = "Waves"             # 1) "Steady current"  2) "Waves"  3) "Tidal current"  4) "Waves + Currents"
#Determineing condition parameters:
U_m, KC, U_c = pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e)
#Determining the scour profile values:
S_eq, W1, W2 = scour_profile(condition, e, D, KC, U_c, U_m)
#plotting the scour profile
Scour_profile_plot(D, W1, W2, S_eq, condition)



#### Waves + Current ####
condition = "Waves + Current"  # 1) "Steady current"  2) "Waves"  3) "Tidal current"  4) "Waves + Current"
#Determining the scour profile values:
S_eq, W1, W2 = scour_profile(condition, e, D, KC, U_c, U_m)
#plotting the scour profile
Scour_profile_plot(D, W1, W2, S_eq, condition)
