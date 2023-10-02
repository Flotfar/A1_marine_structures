
from A1_Input_units import *
from General_constants import *
from Scour_functions import *
from Plots import *


####################################################
#### Input Units
####################################################

#### General units ####
kappa, ks, g, nu, U_f, s, H_s, Fr = general_const(d_50, V, rho_s, rho_w, h, H_s, H)
d = d_50

#### Case specific entities ####
W1_read = 0
W2_read = 0

####################################################
#### Results
####################################################


#### Scour depth developement over time ####

condition = ["Waves", "Steady current"]
# Scour depth at ≈ 100% equibrilum:
S_factor = 0.9999
# Scour depth at 90% and 98% of equilibrium:
S_factor_i = [0.90, 0.98]

#plot collors:
col = ['red', 'blue']
#plot axis labels:
x_lab = "Time [h]"
y_lab = 'Scour depth [m]'
#Labels:
lab = ["Wave scour", "Current scour"]
#other
t_max = 0
s_max = 0


for idx in range(len(condition)):
    if condition[idx] == "Waves":
        #Determineing condition parameters:
        U_m, KC, U_c = pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e)
    elif condition[idx] == "Steady current":
        U, U_cr, U_c = pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g)
        U_m = 0
        KC = 0

    #Determining the scour profile values:
    S_eq, W1, W2 = scour_profile(condition[idx], e, D, KC, U_c, U_m, W1_read, W2_read)
    
    #Determining times at 90% and 98%
    print(condition[idx] + ":")
    for idx2 in range(len(S_factor_i)):
        time = scour_time_scale(U_f, g, s, d, D, h, S_factor_i[idx2], S_eq)
        print("Time to reach " + str(S_factor_i[idx2]*100) + "% scour depth = " + str(round(time, 0)) + " [s]")
        


    time_max = scour_time_scale(U_f, g, s, d, D, h, S_factor, S_eq)
    t = np.linspace(0, time_max, 500)
    
    plt.plot(t/(60*60), scour_depth_time(t, U_f, g, s, d, h, D, S_eq), linestyle='-', color=col[idx], label= lab[idx])

#Grid:
plt.grid(True)
#Legend: 
plt.legend()
#Labels:
plt.xlabel(x_lab)
plt.ylabel(y_lab)

plt.show()
