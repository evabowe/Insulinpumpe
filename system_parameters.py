# Glucose Kinetics
V_G = {"norm_val": 1.88, "diab_val": 1.49}  # dl/kg
k_1 = {"norm_val": 0.065, "diab_val": 0.042}  # min^-1
k_2 = {"norm_val": 0.079, "diab_val": 0.071}  # min^-1

# Insulin Kinetics
V_I = {"norm_val": 0.05, "diab_val": 0.04}  # l/kg
m_1 = {"norm_val": 0.190, "diab_val": 0.379}  # min^-1
m_2 = {"norm_val": 0.484, "diab_val": 0.673}  # min^-1
m_4 = {"norm_val": 0.194, "diab_val": 0.269}  # min^-1
m_5 = {"norm_val": 0.0304, "diab_val": 0.0526}  # min * kg/pmol
m_6 = {"norm_val": 0.6471, "diab_val": 0.8118}  # dimensionless
HE_b = {"norm_val": 0.6, "diab_val": 0.6}  # dimensionless

# Rate of Appearance
k_max = {"norm_val": 0.0558, "diab_val": 0.0465}  # min^-1
k_min = {"norm_val": 0.0080, "diab_val": 0.0076}  # min^-1
k_abs = {"norm_val": 0.057, "diab_val": 0.023}  # min^-1
k_gri = {"norm_val": 0.0558, "diab_val": 0.0465}  # min^-1
f = {"norm_val": 0.90, "diab_val": 0.90}  # dimensionless
a = {"norm_val": 0.00013, "diab_val": 0.00006}  # mg^-1
b = {"norm_val": 0.82, "diab_val": 0.68}  # dimensionless
c = {"norm_val": 0.00236, "diab_val": 0.00023}  # mg^-1
d = {"norm_val": 0.010, "diab_val": 0.09}  # dimensionless

# Endogenous Production
k_p1 = {"norm_val": 2.70, "diab_val": 3.09}  # mg/kg/min
k_p2 = {"norm_val": 0.0021, "diab_val": 0.0007}  # min^-1
k_p3 = {"norm_val": 0.009, "diab_val": 0.005}  # mg/kg/min per pmol/l
k_p4 = {"norm_val": 0.0618, "diab_val": 0.0786}  # mg/kg/min per pmol/kg
k_i = {"norm_val": 0.0079, "diab_val": 0.0066}  # min^-1

# Utilization
F_cns = {"norm_val": 1, "diab_val": 1}  # mg/kg/min
V_m0 = {"norm_val": 2.50, "diab_val": 4.65}  # mg/kg/min
V_mx = {"norm_val": 0.047, "diab_val": 0.034}  # mg/kg/min per pmol/l
K_m0 = {"norm_val": 225.59, "diab_val": 466.21}  # mg/kg

# Secretion
p_2U = {"norm_val": 0.0331, "diab_val": 0.0840}  # min^-1
K = {"norm_val": 2.30, "diab_val": 0.99}  # pmol/kg per (mg/dl)
alpha = {"norm_val": 0.050, "diab_val": 0.013}  # min^-1
beta = {"norm_val": 0.11, "diab_val": 0.05}  # pmol/kg/min per (mg/dl)
gamma = {"norm_val": 0.5, "diab_val": 0.5}  # dimensionless

# Renal Excretion
k_e1 = {"norm_val": 0.0005, "diab_val": 0.0007}  # min^-1
k_e2 = {"norm_val": 339, "diab_val": 269}  # mg/kg

def return_parameters(type):
    if type == "norm":
        return [V_G["norm_val"], k_1["norm_val"], k_2["norm_val"], V_I["norm_val"], m_1["norm_val"], m_2["norm_val"],
                m_4["norm_val"], m_5["norm_val"], m_6["norm_val"], HE_b["norm_val"], k_max["norm_val"],
                k_min["norm_val"], k_abs["norm_val"], k_gri["norm_val"], f["norm_val"], a["norm_val"], b["norm_val"],
                c["norm_val"], d["norm_val"], k_p1["norm_val"], k_p2["norm_val"], k_p3["norm_val"],
                k_p4["norm_val"], k_i["norm_val"], F_cns["norm_val"], V_m0["norm_val"], V_mx["norm_val"],
                K_m0["norm_val"], p_2U["norm_val"], K["norm_val"], alpha["norm_val"], beta["norm_val"],
                gamma["norm_val"], k_e1["norm_val"], k_e2["norm_val"]]
    elif type == "diab":
        return [V_G["diab_val"], k_1["diab_val"], k_2["diab_val"], V_I["diab_val"], m_1["diab_val"], m_2["diab_val"],
                m_4["diab_val"], m_5["diab_val"], m_6["diab_val"], HE_b["diab_val"], k_max["diab_val"],
                k_min["diab_val"], k_abs["diab_val"], k_gri["diab_val"], f["diab_val"], a["diab_val"], b["diab_val"],
                c["diab_val"], d["diab_val"], k_p1["diab_val"], k_p2["diab_val"], k_p3["diab_val"],
                k_p4["diab_val"], k_i["diab_val"], F_cns["diab_val"], V_m0["diab_val"], V_mx["diab_val"],
                K_m0["diab_val"], p_2U["diab_val"], K["diab_val"], alpha["diab_val"], beta["diab_val"],
                gamma["diab_val"], k_e1["diab_val"], k_e2["diab_val"]]