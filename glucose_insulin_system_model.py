# simulated for normal and diabetic type 2
#source of the model:
#Meal Simulation Model of the Glucose-Insulin System Chiara Dalla Man, Robert A. Rizza, and Claudio Cobelli*, Fellow, IEEE

import control as ctrl
import numpy as np
from system_parameters import 

EGP_b = U_b + E_b

dG_p = EGP + Ra - U_ü - E - k1*G_p + k2*G_t     #G_p(0) = G_pb
dG_t = -U_id + k1*G_p - k2*G_t                  #G_t(0) = G_tb
G = G_p / V_G                                  #G(0) = G_b

dI_t = -(m1 +m3)*I_t + m2*I_p +S                #I_t(0) = I_lb
dI_p = -(m2+m4)*I_p + m1*I_l                    #I_p(0) = I_pb
I = I_p / V_I                                   #I(0) = I_b

HE = -m5 * S_t + m6
HE_0 = HE_b
m3 = (HE * m1) / (1 - HE)
m6 = m5 * S_b + HE_b
m3_0 = (HE_b * m1) / (1 - HE_b)
S_b = m3_0 * I_lb + m4 * I_pb
D_b = S_b
m2 = ((S_b / I_pb) - (m4 / (1 - HE_b))) * ((1 - HE_b) / HE_b)
m4 = (2 / 5) * (S_b / I_pb) * (1 - HE_b)

# Endogenous Glucose Production (EGP) equations
EGP = k_p1 - k_p2 * G_p - k_p3 * I_d - k_p4 * I_po
EGP_0 = EGP_b

# Insulin dynamics in the portal vein
dI1 = -k_i * (I1 - I)
dI_d = -k_i * (I_d - I1)
I1_0 = I_b
I_d_0 = I_b

# EGP at basal steady state
k_p1 = EGP_b + k_p2 * G_pb + k_p3 * I_b + k_p4 * I_pob

# Glucose kinetics in the stomach and intestine
Q_sto = Q_sto1 + Q_sto2
dQ_sto1 = -k_gri * Q_sto1 + D*d
dQ_sto2 = -k_empt * Q_sto2 + k_gri*Q_sto1
dQ_gut = -k_abs*Q_gut +k_epmt*Q_sto2
Ra = (f*k_abs*Q_gut)/BW 

Q_sto0 = 0
Q_sto10 = 0
Q_sto20 = 0
Q_gut0 = 0
Ra0 = 0

EGP_b = F_cns + (V_m0 * G_tb)/(K_m0+G_tb)
U_b = EGP_b 

V_m0 = ((EGP_b - F_cns)*(K_m0 + G_tb))/G_tb

U_ii = F_cns

U_id = (V_m * G_t)/(K_m + G_t)
V_m = V_m0 + V_mx * X
K_m = K_m0 + K_mx * X

dX = -p_2U * X + p_2U * (I - I_b)  #! check
X0 = 0

U = U_ii + U_id

G_tb = (F_cns - EGP_b + k1*G_pb)/ k2

V_m0 = ((EGP_b - F_cns)* (K_m0 + G_tb))/G_tb
S = gamma * I_po
dI_po = - gamma* I_po + S_po
I_po0 = I_pob

S_po = Y + K * dG + S_b if dG > 0 else Y + S_b          #! check

E = k_e1 * (G_p - k_e2) if G_p > k_e2 else 0

dY = -alpha * (Y- beta*(G-h)) if beta*(G-h) >= -S_b else -alpha*Y - alpha * S_b
Y0 = 0


