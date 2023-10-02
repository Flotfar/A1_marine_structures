
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

#plot collors:
col = ['red', 'blue']
#plot axis labels:
x_lab = "Time [h]"
y_lab = 'Scour depth [m]'

#other
t_max = 0
s_max = 0


for idx in range(2):

    ##### Waves backfilling input #####
    U_m, KC, U_c = pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e)
    S_eq_wave, W1, W2 = scour_profile(condition[0], e, D, KC, U_c, U_m, W1_read, W2_read)

    ##### Current input #####
    U, U_cr, U_c = pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g)
    U_m = 0
    KC = 0
    S_eq_current, W1, W2 = scour_profile(condition[1], e, D, KC, U_c, U_m, W1_read, W2_read)
    
    time_start = 14 * idx
    time_end = 14 * (1+idx)
    
    t_time = np.linspace(time_start, time_end, 500)
    t = np.linspace(0, (14 * 3600), 500)

    if idx == 0:
        S_eq = S_eq_current
        y_val = scour_depth_time(t, U_f, g, s, d, h, D, S_eq)
    else:
        S_eq = S_eq_wave
        S_i = scour_depth_time((14 * 3600), U_f, g, s, d, h, D, S_eq_current)
        y_val =  backfilling(U_f, g, s, d, S_i, S_eq, h, D, t)
    
    plt.plot(t_time, y_val, linestyle='-', color=col[idx])


#Grid:
plt.grid(True)
#Legend: 
plt.legend()
#Labels:
plt.xlabel(x_lab)
plt.ylabel(y_lab)

plt.show()




