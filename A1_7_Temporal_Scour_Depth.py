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

#Finding shielding parameter from current+wave velocity ucw:
    
    
U_dfwc = Ucwd/(8.6+2.5*np.log(D/(2*kb)))
theta_dcw = (U_dfwc**2)/(g*(s-1)*d50)


#Finding gamma from slide 53:

gamma = 0.31 #Since KC<2.

#Finding alpha from slide 53:
Alpha = 1-0.9*np.tanh(3000*Ucwd**20-0.6**20)   
 
#Current to waves plus current
Tstarb = gamma*Alpha*theta_dcw**(-3/2)

# (2.44) on slide 44
T = (D**2/(g*(s-1)*d50**3)**(0.5))*Tstarb

#Plotting scour with backfilling from lecture 1 slide 59

def Seq(t):
    return S_dwceq+(S_swceq-S_dwceq)*np.exp(-t/Tstarb)

t = np.linspace(0,300,1000)

plt.plot(t,Seq(t),color='red')
plt.ylim(0,6)
plt.ylabel('Scour depth, S (m)')
plt.xlabel('Time, t (min)')
plt.savefig('plot.png', dpi=300)