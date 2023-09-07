import numpy as np


##### Pipeline - Scour for a steady current: #####
def pipeline_scour_steady_current(D, e, U_f, kappa, ks, s, n, g):

    z = D - e
    U = U_f/kappa * np.log((30*z)/ks)
    U_cr = np.sqrt(0.025*np.exp(9*(0.05/D)**(0.5))*g*D*(1-n)*(s-1))

    #Evalutaing the onset criterias
    if U >= U_cr:
        print("\n" + str(round(U, 3)) + ">=" + str(round(U_cr, 3)))
        print("Scour can occour around the pipeline, from steady current \n")
    else:
        print("\n" + str(round(U, 3)) + "<" + str(round(U_cr, 3)))
        print("Scour cannot occour around the pipeline, from steady current \n")

    return U, U_cr


##### Pipeline - Scour from waves: #####
def pipeline_scour_waves(T_p, T, D, g, H_s, h, H, L, t, x, U_f, ks, kappa, e):

    z = D
    U = U_f/kappa * np.log((30*z)/ks)

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

    return U_m, KC

def pipeline_scour_waves_Ucr(reading, g, D, n, s):
    U_mcr = np.sqrt(g*D*(1-n)*(s-1)*reading)
    print("U_mcr = " + str(round(U_mcr, 3)) + "\n")

    return U_mcr


#### Migration of span shoulders ####
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
    S_h_div = D*(4)
    # S_h_div = D*(15+4)

    time_Vh1 = S_h/V_h1
    time_Vh2 = S_h/V_h2

    time_div_Vh1 = (S_h_div/V_h1)
    time_div_Vh2 = (S_h_div/V_h2)
    # time_div_Vh1 = (S_h_div/V_h1)
    # time_div_Vh2 = (S_h_div/V_h2)

    time = time_Vh2 - time_Vh1
    time_div = time_div_Vh2 - time_div_Vh1

    print("Migration time from primary to secondary: " + str(round(time,3)) + " ± " + str(round(time_div, 3))+ " [s] \n")
    #print("Migration time:  " + str(round(time_Vh1, 3) + "±" + str(time_div_Vh1, 3)))

    return V_h1, V_h2, time, time_div
