import numpy as np

#general expresions
def general_const(d_50, V, rho_s, rho_w, h, H_s, H):

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
    return kappa, ks, g, nu, U_f, s, H_s, Fr