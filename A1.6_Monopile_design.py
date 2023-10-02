
import numpy as np
from scipy.optimize import fsolve

from A1_Mono_Input_units import *
from General_constants import *
from Scour_functions import *


####################################################
#### Input Units
####################################################

#### General units ####
V, H_s, H = 0, 0, 0
kappa, ks, g, nu, U_f, s, H_s, Fr = general_const(d_50, V, rho_s, rho_w, h, H_s, H)

####################################################
#### Calculations and results
####################################################


####################################################
#### a)
####################################################
#### Summer conditions: ####

#Resulting in U/Ucr = 1.77 and by reading from figure 3.25 we get:
KI_s = 1.6
#For boundary layer thickness we assume delta=h, resulting in delta/D = 3.33
#From 3.26 we then get:
Ks_delta = 0.95
# For grain size D/d_50 = 45000 >> 30 meaning according to fig 3.27:
Ks_d = 1
#The cross section shape is assumed ciruclar, resulting in:
Ks_s = 1.0
#Allignment is 0 degrees, meaning:
Ks_alpha = 1

U_sc, S_seq, U_s, Rey_s = eq_scour_current(V, d_50, d_84, g, s, h, nu, KI_s, Ks_delta, Ks_d, Ks_s, Ks_alpha, D)

print("U_sc = " + str(round(U_sc,3)) +"\nU_s = " + str(round(U_s,3)) + "\ngrain Reynolds number = " + str(round(Rey_s,3)))
print("The equilibrium scour depth for the summer conditions is: " + str(round(S_seq,3)))



#### Design conditions: ####

#Resulting in U/Ucr = 2.78 and by reading from figure 3.25 we get:
KI_d = 1.4
#For boundary layer thickness we assume delta=h, resulting in delta/D = 3.33
#From 3.26 we then get:
Kd_delta = 0.95
# For grain size D/d_50 = 45000 >> 30 meaning according to fig 3.27:
Kd_d = 1
#The cross section shape is assumed ciruclar, resulting in:
Kd_s = 1.0
#Allignment is 0 degrees, meaning:
Kd_alpha = 1

U_d, Rey_d, U_dc, S_deq = eq_scour_current(V_d, d_50, d_84, g, s, h, nu, KI_d, Kd_delta, Kd_d, Kd_s, Kd_alpha, D)

print("\nU_sc = " + str(round(U_dc,3)) +"\nU_s = " + str(round(U_d,3)) + "\ngrain Reynolds number = " + str(round(Rey_d,3)))
print("The equilibrium scour depth for the design conditions is: " + str(round(S_deq,3)))



####################################################
#### b)
####################################################

#### Summer conditions: ####

#Calculating the free-stream velocity magnitude for irregular waves (ASSUMING CLEAR WATER)
S_swceq, Ucws = eq_scour_wavecurrent(T_ps, H_ss, g, h, U_s, D, S_seq)

print("\nflow ratio = " + str(round(Ucws,3)))
print("The equilibrium scour depth for the summer conditions with combined waves and currents is:")
print(round(S_swceq,3))


#### Design Conditions: ####
    
S_dwceq, Ucwd = eq_scour_wavecurrent(T_pd, H_sd, g, h, U_d, D, S_deq)

print("\nflow ratio = " + str(round(Ucwd,3)))
print("The equilibrium scour depth for the design conditions with combined waves and currents is:")
print(round(S_dwceq,3))

