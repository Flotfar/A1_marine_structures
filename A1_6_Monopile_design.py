import numpy as np

from A6_Input_units import *
from General_constants import *
from Scour_functions import *
from scipy.optimize import fsolve

####################################################
#### Input Units
####################################################

# General units (Summer)
D=9 #m
h=30 #m
d50=0.2*10**(-3) #m
sigmag=2.3 #-
rhos=2650 #kg/m3
rho=1026 #kg/m3
kappa = 0.4
kb = 2.5*d50
g = 9.81
s = rhos/rho
nu = 10**(-6)

Vs = 0.7 #m/s
Hss = 3.5 #m
Tps = 8.0 #s

# General units (Design)
Vd = 1.1 #m/s
Hsd = 8.8 #m
Tpd = 14.0 #s

####################################################
#### Calculations and results
####################################################

###########
# a)
###########

#summer:

Us = Vs

#Finding Grain reynold's number:
Rey_s = Us*d50/nu

#From figure 1.2 the critical shield parameter can be found from the grain
# reynold's number:
theta_scr = 4.5*10**(-2)

#From (1.3) the critical friction velocity is found by solving for Uf:
U_sfcr = np.sqrt(theta_scr*g*d50*s-theta_scr*g*d50)

# From the formula for mean depth averaged current we solve for V, getting the critical flow velocity Uc:

U_sc = (6.0+2.5*np.log(h/kb))*U_sfcr

#Resulting in U/Ucr = 1.77 and by reading from figure 3.25 we get:
KIs = 1.6

#For boundary layer thickness we assume delta=h, resulting in delta/D = 3.33
#From 3.26 we then get:
Ksdelta = 0.95

# For grain size D/d50 = 45000 >> 30 meaning according to fig 3.27:
Ksd = 1

#The cross section shape is assumed ciruclar, resulting in:
Kss = 1.0

#Allignment is 0 degrees, meaning:
Ksalpha = 1

#The equilibrium scour depth is then:
S_seq = KIs*Ksdelta*Ksd*Kss*Ksalpha*D

print("The equilibrium scour depth for the summer conditions is:")
print(S_seq)

# Design

Ud = Vd

#Finding Grain reynold's number:
Rey_d = Ud*d50/nu

#From figure 1.2 the critical shield parameter can be found from the grain
# reynold's number:
theta_dcr = 4.5*10**(-2)

#From (1.3) the critical friction velocity is found by solving for Uf:
U_dfcr = np.sqrt(theta_dcr*g*d50*s-theta_dcr*g*d50)

# From the formula for mean depth averaged current we solve for V, getting the critical flow velocity Uc:

U_dc = (6.0+2.5*np.log(h/kb))*U_dfcr

#Resulting in U/Ucr = 2.78 and by reading from figure 3.25 we get:
KId = 1.4

#For boundary layer thickness we assume delta=h, resulting in delta/D = 3.33
#From 3.26 we then get:
Kddelta = 0.95

# For grain size D/d50 = 45000 >> 30 meaning according to fig 3.27:
Kdd = 1

#The cross section shape is assumed ciruclar, resulting in:
Kds = 1.0

#Allignment is 0 degrees, meaning:
Kdalpha = 1

#The equilibrium scour depth is then:
S_deq = KId*Kddelta*Kdd*Kds*Kdalpha*D

print("The equilibrium scour depth for the design conditions is:")
print(S_deq)

###########
# b)
###########

#Summer:

#Calculating the free-stream velocity magnitude for irregular waves (ASSUMING CLEAR WATER)
Tzs = Tps/1.3
Ums = (Hss/(2*np.sqrt(2)))*np.sqrt(g/h)*np.exp(-((3.65/Tzs)*np.sqrt(h/g))**2.1)

Ucws = Us/(Us+Ums)

# Keulegan-Carpenter number:
KCs = Ums*Tps/D

#Calculating A and B:

#Determining L first:

#Deep water wavelength:
L0s = g*Tps**2/(2*np.pi)

#Corresponding wave number:
k0s = 2 * np.pi/L0s

#solve for actual wave number k:
    
def equation(k):
    return g * k * np.tanh(k * h) - (2 * np.pi / Tps)**2

ks = fsolve(equation, k0s)

# Calculating the wave length
Ls = float(2 * np.pi / ks)

As = 0.03 + 8*Ucws**(1/max(KCs,0.5)+5)
Bs = (6-5.8*np.tanh(200*(D/Ls)**1.9))*np.exp(-4.7*Ucws)

#Eq scour depth for summer conditions for combined wave and current:
    
S_swceq = ((S_seq/D)*(1-np.exp(-As*(max(KCs,0.5)-Bs))))*D

print("The equilibrium scour depth for the summer conditions with combined waves and currents is:")
print(S_swceq)

#Design:
    
#Calculating the free-stream velocity magnitude for irregular waves (ASSUMING CLEAR WATER)
Tzd = Tpd/1.3
Umd = (Hsd/(2*np.sqrt(2)))*np.sqrt(g/h)*np.exp(-((3.65/Tzd)*np.sqrt(h/g))**2.1)

Ucwd = Ud/(Ud+Umd)

# Keulegan-Carpenter number:
KCd = Umd*Tpd/D

#Calculating A and B:

#Determining L first:

#Deep water wavelength:
L0d = g*Tpd**2/(2*np.pi)

#Corresponding wave number:
k0d = 2 * np.pi/L0d

#solve for actual wave number k:
    
def equation(k):
    return g * k * np.tanh(k * h) - (2 * np.pi / Tpd)**2

kd = fsolve(equation, k0d)

# Calculating the wave length
Ld = float(2 * np.pi / kd)

Ad = 0.03 + 8*Ucwd**(1/max(KCd,0.5)+5)
Bd = (6-5.8*np.tanh(200*(D/Ld)**1.9))*np.exp(-4.7*Ucwd)

#Eq scour depth for summer conditions for combined wave and current:
    
S_dwceq = ((S_deq/D)*(1-np.exp(-Ad*(max(KCd,0.5)-Bd))))*D

print("The equilibrium scour depth for the design conditions with combined waves and currents is:")
print(S_dwceq)