####################################################
#### Input Units
####################################################

#Monopile foundation
D = 9.0  #[m] foundation diameter

#Bed conditions
h = 30   #[m] water depth around monopile
d_50 = 0.2  #[mm] grain size
d_84 = 2.3*d_50  #[mm] grain size
sigma_g = d_84/d_50
rho_s = 2650   #[kg/m3] sand density
rho_w = 1026   #[kg/m3] sea-water density

#Summer and design conditions
V_s = 0.7  #[m/s]
V_d = 1.1  #[m/s]
H_ss = 3.5 #[m]
H_sd = 8.8 #[m]
T_ps = 8.0 #[s]
T_pd = 14.0 #[s]
