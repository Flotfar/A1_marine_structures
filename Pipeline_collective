import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


####################################################
#### Input Units
####################################################

d_50 = 0.2 * 10**(-3)  #[m] input in mm
rho_s = 2650        #Sand Density [kg/m3] 
rho_w = 1026        #Water Density [kg/m3]
n = 0.43            #Porosity [-]
V = 1.1             #Dominant velocity [m/s]
h = 15              #Water depth [m]
#pipeline input:
D = 0.8             #Diameter of the pipeline [m]
e = 4 *10**(-2)     #Initial Burial depth [m] input in cm
#regular waves:
T_p = 9.0           #Peak wave period (irregular waves) [s]
H_s = 4.0           #Significant wave height [m]
# regular waves:
T = 0               #Wave period (regular waves) [s]
H = 0               #Wave height [m]
L = 0               #Wave lenght [m]
t = 0
x = 0


####################################################
#### General constants
####################################################

kappa = 0.4             #Von karman constant [-]
ks = 2.5*d_50           #Nikuradse’s bed roughness [m]
g = 9.81                #Gravity in denmark [m/s2]
nu = 1.007 * 10**(-6)               #Viscosity at 20 deg [m2/s]
U_f = V/(6 + 2.5 * np.log(h/ks))    #Friction Velocity - steady current []
Fr = V/(np.sqrt(g*h))   #Froude number


#Relative density [-]
if rho_s != 0 and rho_w != 0:       
    s = rho_s/rho_w
else:
    s = 2.65

#Significant waveheight:
if H_s == 0:
    H_s = H*(1/3)




####################################################
#### Functions
####################################################


##### Pipeline - Scour for a steady current: #####
# Return:: U, U_cr, U_c
def pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g):

    z = D - e
    z2 = D/2
    U = U_f/kappa * np.log((30*z)/ks)
    U_c = U_f/kappa * np.log((30*z2)/ks)
    U_cr = np.sqrt(0.025*np.exp(9*(0.05/D)**(0.5))*g*D*(1-n)*(s-1))
    
    #Evalutaing the onset criterias
    if U >= U_cr:
        print("\n" + str(round(U, 3)) + ">=" + str(round(U_cr, 3)))
        print("Scour can occour around the pipeline, from steady current \n")
    else:
        print("\n" + str(round(U, 3)) + "<" + str(round(U_cr, 3)))
        print("Scour cannot occour around the pipeline, from steady current \n")

    return U, U_cr, U_c


##### Pipeline - Scour from waves: #####
# Return: U_m, KC, U_c
def pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e):

    z = D/2
    U_c = U_f/kappa * np.log((30*z)/ks)

    #Irregular waves:
    if T_p != 0:
        T_z = T_p/1.3
        U_m = H_s/(2*np.sqrt(2))*np.sqrt(g/h)*np.exp(-((3.65/T_z)*np.sqrt(h/g))**2.1) #Laminar boundary layer
        KC = (U_m * T_p)/(D)                        #Keulegan-Carpenter number (Wave stroke to diameter ratio)

        print("U_m = " + str(str(round(U_m, 3))) + " | " + "KC = " +str(str(round(KC, 3))) + " | " + "e/D = " +str(str(round(e/D, 3))))
        print("U_cr must be found from fig.2.8 Sumer and fredsoe 2002 \n")

    #Regular waves:
    if T != 0:
        
        Omega = (2*np.pi)/T         #Wave frequency
        k = (2 * np.pi)/(L)         #Wave number
        U_m = (np.pi * H)/(T) * (1)/(np.sinh(k * h)) #Laminar boundary layer
        U_0 = U_m * np.cos(Omega * t - k * x)        #Near bed orbital velocity (dont know  meaning of t and x)
        KC = (U_m * T)/(D)          #Keulegan-Carpenter number (Wave stroke to diameter ratio)

        print("U_m = " + str(str(round(U_m, 3))) + " | " + "KC = " +str(str(round(KC, 3))) + " | " + "e/D = " +str(str(round(e/D, 3))))
        print("U_cr must be found from fig.2.8 Sumer and fredsoe 2002 \n")

    return U_m, KC, U_c

# Return: U_mcr
def pipeline_scour_waves_Ucr(reading, g, D, n, s):
    U_mcr = np.sqrt(g*D*(1-n)*(s-1)*reading)
    print("U_mcr = " + str(round(U_mcr, 3)) + "\n")

    return U_mcr


#### Migration of span shoulders ####
# Return: V_h1, V_h2, time
def migration_span_shoulders(regime, U_f, nu, g, s, d, e, D):

    Re_g = (U_f * d) / (nu)             # grain Reynolds number
    theta = (U_f**2) / (g * (s-1) * d)  # Shields parameter
    theta_c = 0.165*(Re_g + 0.6)**(-0.8) + 0.045 * np.exp(-40*Re_g**(-1.3)) #critical shields parameter
    
    #### Migration speed ####
    if regime == "live_bed":
        
        V_h11 = 3*(1+200*(theta-theta_c)**(3/2))*np.exp(-3.2*(e/D))
        V_h21 = 3*(1+200*(theta-theta_c)**(3/2))*np.exp(-6*(e/D))

        V_h1 = V_h11 * np.sqrt(g*(s-1)*d)*(d/D)
        V_h2 = V_h21 * np.sqrt(g*(s-1)*d)*(d/D)
    
    elif regime == "clear_water":

        V_h11 = 3*(np.tanh(4.2*((theta/theta_c)-0.12)**(3.5)))*np.exp(-3.2*(e/D))
        V_h1 = V_h11 * np.sqrt(g*(s-1)*d)*(d/D)
        V_h2 = 0

    print("Primary migration speed = " + str(round(V_h1, 3)) + " [m/s] \n" + "Secondary migration speed = " + str(round(V_h2, 3))+ " [m/s]")

    #### Migration time ####
    S_h = D*15
    S_v = D*0.6
    S_h_div = D*(4)
    S_v_div = D*(0.2)
    
    time_h = S_h/V_h1
    time_v = S_v/V_h1
    
    time_div_h = (S_h_div/V_h1)
    time_div_v = (S_v_div/V_h1)

    if time_h > time_v:
        time = time_h
        time_div = time_div_h
        print("S_h is dominant")
    else:
        time = time_v
        time_div = time_div_v
        print("S_v is dominant")
   
    print("Migration time from primary to secondary: " + str(round(time,3)) + " ± " + str(round(time_div, 3))+ " [s] \n")

    return V_h1, V_h2, time


#### Scour profile values ####
# Return: S_eq, W1, W2
def scour_profile(condition, e, D, KC, U_c, U_m, W1_read, W2_read):
    if condition == "Steady current":
        #Equlibrium Scour depth "S_eq"
        rel = e/D
        if rel >= -0.25 and rel <= 1.2:
            S_eq = D * 0.625 * np.exp(-0.6 * e/D)
        else:
            print("Steady current condition e/D are not met")
        #Scour Width "W"
        W1 = 2*D
        W2 = 4*D
    
    elif condition == "Tidal current":
        #Equlibrium Scour depth "S_eq"
        S_eq = 2.6*D
        #Scour Width "W"
        W1 = 0.35 * KC**(0.65) * D
        W2 = 0.35 * KC**(0.65) * D


    elif condition == "Waves":
        #Equlibrium Scour depth "S_eq"
        rel = e/D
        if rel >= 0 and rel <= 2:
            S_eq = 0.1 * np.sqrt(KC) * np.exp(-0.6 * e/D) * D
        else:
            print("Wave condition e/D are not met")
        #Scour Width "W"
        W1 = 0.35 * KC**(0.65) * D
        W2 = 0.35 * KC**(0.65) * D


    elif condition == "Waves + Current":

        rel2 = U_c / (U_c + U_m)

        #Equlibrium Scour depth "S_eq"
        rel1 = e/D
        if rel1 >= -0.25 and rel1 <= 1.2:
            S_cur = D * 0.625 * np.exp(-0.6 * e/D)
        else:
            print("Steady current condition e/D are not met")
        
        if 0 < rel2 and rel2 <= 0.7:

            if rel2 > 0 and rel2 <= 0.4:
                a = 0.557 - 0.912 * (rel2 - 0.25)**2
                b = -1.14 + 2.24 * (rel2 - 0.25)**2
                F = 5/3 * (KC)**a * np.exp(2.3*b)
            
            elif 0.4 < rel2 and rel2 <= 0.7:
                a = -2.14 * rel2 + 1.46
                b = 3.3 * rel2 - 2.5
                F = 5/3 * (KC)**a * np.exp(2.3*b)
        
        elif 0.7 < rel2 and rel2 <= 1:
            F = 1
        else: 
            print("Error with Uc/(Uc + Um) condition")
        
        S_eq = S_cur * F
    
        #Scour Width "W"
        W1 = W1_read * D
        W2 = W2_read * D


    return S_eq, W1, W2


#### Scour time scale ####
# Return: t
def scour_time_scale(U_f, g, s, d, D, h, S_factor, S_eq):

    theta = (U_f**2) / (g * (s-1) * d)  # Shields parameter
    T_nondim = (1/375) * (h/D)**(0.75) * theta**(-1.5)  # Non-dimensional time scale
    T = (D**2)/(np.sqrt(g*(s-1)*d**3)) * T_nondim  # Time scale
    S_t = S_eq * S_factor  # The scour depth at given factor of the eq-scour depth
    
    time = -np.log((S_eq - S_t)/(S_eq)) * T  # time for scouring to develope

    return time

def scour_depth_time(t, U_f, g, s, d, h, D, S_eq):

    theta = (U_f**2) / (g * (s-1) * d)  # Shields parameter
    T_nondim = (1/375) * (h/D)**(0.75) * theta**(-1.5)  # Non-dimensional time scale
    T = (D**2)/(np.sqrt(g*(s-1)*d**3)) * T_nondim  # Time scale

    S_time = S_eq * (1-np.exp(-t/T))  # Scour depth to the time t

    return S_time

def backfilling(U_f, g, s, d, S_i, S_eq, h, D, t):

    theta = (U_f**2) / (g * (s-1) * d)
    T_nondim = (1/375) * (h/D)**(0.75) * theta**(-1.5)  # Non-dimensional time scale
    T = (D**2)/(np.sqrt(g*(s-1)*d**3)) * T_nondim  # Time scale
    Tb = T * 0.2*theta**(-5/3)

    S_time = S_eq + (S_i - S_eq)*np.exp(-t/Tb)

    return S_time


####################################################
#### Plot functions
####################################################

#### Scour profile plot
def Scour_profile_plot(D, W1, W2, S_eq, condition):
    # Coordinates of three points (W1, W2, and S)
    p1 = (-W1, 0)
    p2 = (W2, 0)
    p3 = (0, -S_eq)

    # Calculate the axis limits based on W1 and W2
    factor = 0.2
    x_min = min(-W1 - W1*factor , -2)
    x_max = max(W2 + W2*factor, 4)
    y_min = min(-S_eq - S_eq*factor, -x_max/2, -2)
    y_max = max(x_max/2, 2)

    # Create a fixed-size figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))  # Adjust the size as needed

    # Create a circle patch
    circle = Circle((0, D / 2), radius=D / 2, fill=False, color='black')  # (0, D/2) is the center of the circle

    # Add the circle patch to the axis
    ax.add_patch(circle)

    # Plot the three input points
    ax.plot(*p1, marker='o', markersize=3, color='black', label='W1')
    ax.plot(*p2, marker='o', markersize=3, color='black', label='W2')
    ax.plot(*p3, marker='o', markersize=3, color='black', label='S')

    # Set the aspect ratio to 'equal' to ensure the circle looks circular
    ax.set_aspect('equal', adjustable='box')

    # Set axis limits to include the circle
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Add a grid
    ax.set_axisbelow(True)
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7, zorder=0)

    # Plot a solid line for y=0
    plt.axhline(0, color='black', linestyle='dashed', linewidth=1)

    # Optional: Add labels or customize the plot further
    ax.set_xlabel('Lenght [m]')
    ax.set_ylabel('Height [m]')
    ax.set_title('Scour profile in ' + condition)

    # Display the plot
    plt.show()




####################################################
#### Question 1
####################################################

print("\n######## Results: ########")
#### Scour for a steady current ####
U, U_cr = pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g)

#### Scour from waves: ####
U_m, KC, U_c = pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e)

print("##########################\n")

# when gaining the result input the reading (input a pop up menu)
reading = 0.05
U_mcr = pipeline_scour_waves_Ucr(reading, g, D, n, s)



####################################################
#### Question 2
####################################################

#### Case specific units ####
regime = "live_bed"    #either "live_bed" or "clear_water"

print("\n######## Results: ########")

#### Migration speeds and time from primary to secondary ####
d = d_50
V_h1, V_h2, time = migration_span_shoulders(regime, U_f, nu, g, s, d, e, D)



####################################################
#### Question 3
####################################################

#### Case specific entities ####
W1_read = 0
W2_read = 0
U_m = 0
KC = 0

#### Current alone ####
condition = "Steady current"    # 1) "Steady current"  2) "Waves"  3) "Tidal current"  4) "Waves + Currents"
#Determineing condition parameters:
U, U_cr, U_c = pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g)
#Determining the scour profile values:
S_eq, W1, W2 = scour_profile(condition, e, D, KC, U_c, U_m, W1_read, W2_read)
#plotting the scour profile
Scour_profile_plot(D, W1, W2, S_eq, condition)



#### Tidal current ####
condition = "Tidal current"    # 1) "Steady current"  2) "Waves"  3) "Tidal current"  4) "Waves + Currents"
#Determining the scour profile values:
KC = 990
S_eq, W1, W2 = scour_profile(condition, e, D, KC, U_c, U_m, W1_read, W2_read)
#plotting the scour profile
Scour_profile_plot(D, W1, W2, S_eq, condition)




####################################################
#### Question 4
####################################################

#### Case specific entities ####
W1_read = 0
W2_read = 0
d = d_50

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



####################################################
#### Question 5
####################################################
#### Case specific entities ####
W1_read = 0
W2_read = 0

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
