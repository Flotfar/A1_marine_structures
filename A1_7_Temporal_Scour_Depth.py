import numpy as np

from A6_Input_units import *
from A1_6_Monopile_design import *
from General_constants import *
from Scour_functions import *
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

####################################################
#### Calculations and results
####################################################

#Since summer depth is deeper than design depth, backfilling occurs in the shift.

#Slide 53

#Finding shield parameter for current:

U_fc = Ud/(6+2.5*np.log(h/(2.5*d50)))
theta_c = (U_fc**2)/(g*(s-1)*d50)


Tz = Tpd/1.3
Um=Hsd/(2*np.sqrt(2))*np.sqrt(g/h)*np.exp(-(3.65/Tz*np.sqrt(h/g))**2.1)

a = Um*Tpd/(2*np.pi)
fw_rough = np.exp(5.5*(a/(2.5*d50))**(-0.16)-6.7)
Rew = a*Um/nu

fw_smooth = 0.035*Rew**(-0.16)

fw_laminar = 2/np.sqrt(Rew)

fw = max(fw_rough,fw_smooth,fw_laminar)

U_fw = np.sqrt(fw/2)*Umd

theta_w = (U_fw**2)/(g*(s-1)*d50)

theta_m = theta_c*(1+1.2*(theta_w/(theta_c+theta_w))**3.2)

theta_cw = theta_m + theta_w

Uc = (U_fc/0.4)*np.log(30*(D/2)/(2.5*d50))
m = Uc/(Uc+Umd)

#Finding gamma from slide 53:
KC = (Um*Tpd)/D

gamma = 0.31 #Since KC<2.

#Finding alpha from slide 53:
Alpha = 1-0.9*np.tanh(3000*m**20-0.6**20)   
 
#Current to waves plus current
Tstarb = gamma*Alpha*theta_cw**(-3/2)

# (2.44) on slide 44
T = (D**2/(g*(s-1)*d50**3)**(0.5))*Tstarb

#Plotting scour with backfilling from lecture 1 slide 59

def Seq(t):
    return S_dwceq+(S_swceq-S_dwceq)*np.exp(-t/T)

t = np.linspace(0,9*10**6,1000)

plt.plot(t,Seq(t),color='red')
plt.ylim(0,6)
plt.ylabel('Scour depth, S (m)')
plt.xlabel('Time, t (s)')
plt.savefig('plot.png', dpi=300)