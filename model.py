from system_parameters import return_parameters
import numpy as np

def G_p(iteration):
    G_p0 = G_pb
    G_p = G_p0
    for i in range(iteration):
        dG_p = EGP(i) + Ra(i) - U_ii() - E(i) - k_1*G_p + k_2*G_t(i)
        G_p = G_p + dG_p
    return G_p
def U_ii():
    return F_cns
def G_t(iteration):
    G_t0 = G_tb
    G_t = G_t0
    for i in range(iteration):
        dG_t = -U_id(i) + k_1*G_p(i) - k_2*G_t
        G_t = G_t + dG_t
    return G_t
def G(iteration):
    G0 = G_b
    G = G0
    for i in range(iteration):
        G = G_p(i) / V_G
    return G

def U_id(iteration):
    i = iteration
    U_id = (V_m(i) * G_t(i))/(K_m(i) + G_t(i))
    return U_id

def K_m(iteration):
    i = iteration
    K_m = K_m0  # + K_mx * X(i)
    return K_m

def V_m(iteration):
    i = iteration
    V_m = V_m0 + V_mx * X(i)
    return V_m

def X(iteration):
    X0 = 0
    X = X0
    for i in range(iteration):
        dx = -p_2U * X + p_2U * (I(i) - I_b)
        X = X + dx
    return X


def I_l(iteration):
    I_l0 = I_lb
    I_l = I_l0
    for i in range(iteration):
        dI_l = -(m_1 +m_3(i))*I_l + m_2*I_p(i) + S(i)
        I_l = I_l + dI_l
    return I_l
def S(iteration):
    i = iteration
    S = gamma * I_po(i) + S_po(i)
    return S
def Ra(iteration):
    BW = 75
    i = iteration
    Ra = (f*k_abs*Q_gut(i))/BW 
    return Ra

def E(iteration):
    i = iteration
    E = k_e1 * (G_p(i) - k_e2) if G_p(i) > k_e2 else 0
    return E

def Q_gut(iteration):
    Q_gut0 = 0
    Q_gut = Q_gut0
    for i in range(iteration):
        dQ_gut = -k_abs*Q_gut +k_empt(i)*Q_sto2(i)
        Q_gut = Q_gut + dQ_gut
    return Q_gut

def Q_sto2(iteration):
    Q_sto20 = 0
    Q_sto2 = Q_sto20
    for i in range(iteration):
        dQ_sto2 = -k_empt(i) * Q_sto2 + k_gri*Q_sto1(i)
        Q_sto2 = Q_sto2 + dQ_sto2
    return Q_sto2
def Q_sto(iteration):
    Q_sto = Q_sto1(iteration) + Q_sto2(iteration)
    return Q_sto


def k_empt(iteration):
    i = iteration
    k_empt = k_min + ((k_max - k_min) / 2) * (np.tanh(alpha * (Q_sto(i) - b*D())) - np.tanh(beta * (Q_sto(i) - c*D())) +2)
    return k_empt

def Q_sto1(iteration):
    Q_sto10 = 0
    Q_sto1 = Q_sto10
    for i in range(iteration):
        dQ_sto1 = -k_gri * Q_sto1 + D()*d
        Q_sto1 = Q_sto1 + dQ_sto1
    return Q_sto1

def D():            #! Ingested glucose
    return 70


def I_p(iteration):
    I_p0 = I_pb
    I_p = I_p0
    for i in range(iteration):
        dI_p = -(m_2+m_4)*I_p + m_1*I_l(i)
        I_p = I_p + dI_p
    return I_p

def I(iteration):
    I0 = I_b
    I = I0
    for i in range(iteration):
        I = I_p(i) / V_I
    return I
def I1(iteration):
    I1_0 = I_b
    I1 = I1_0
    for i in range(iteration):
        dI1 = -k_i * (I1 - I(i))
        I1 = I1 + dI1
    return I1

def m_3(iteration):
    i = iteration
    m_3 = (HE(i) * m_1) / (1 - HE(i))
    return m_3

def HE(iteration):
    i = iteration
    HE = -m_5 * S(i) + m_6
    return HE

def I_d(iteration):
    I_d_0 = I_b
    I_d = I_d_0
    for i in range(iteration):
        dI_d = -k_i * (I_d - I1(i))
        I_d = I_d + dI_d
    return I_d
def I_po(iteration):
    I_po0 = I_pob
    I_po = I_po0
    for i in range(iteration):
        dI_po = - gamma* I_po + S_po(i)
        I_po = I_po + dI_po
    return I_po

def S_po(iteration):
    i = iteration
    dG = G_p(i) - G_p(i-1)
    S_po = Y(i) + K * dG + S_b() if dG > 0 else Y(i) + S_b()
    return S_po

def Y(iteration):
    Y = 0
    h = G(0)
    for i in range(iteration):
        dY = -alpha * (Y- beta*(G(i)-h)) if beta*(G(i)-h) >= -S_b() else -alpha*Y - alpha * S_b()
        Y = Y + dY
    return Y

def S_b():
    return m3_0() * I_lb + m_4 * I_pb

def m3_0():
    return (HE_b * m_1) / (1 - HE_b)

def EGP(iteration):
    EGP_b = F_cns + (V_m0 * G_tb)/(K_m0+G_tb)
    EGP_0 = EGP_b
    EGP = EGP_0
    for i in range(iteration):
        EGP = k_p1 - k_p2 * G_p(i) - k_p3 * I_d(i) - k_p4 * I_po(i)
    return EGP

if __name__ == "__main__":
    global V_G, k_1, k_2, V_I, m_1, m_2, m_4, m_5, m_6, HE_b, k_max, k_min, k_abs, k_gri, f, a, b, c, d, k_p1, k_p2, k_p3, k_p4, k_i, F_cns, V_m0, V_mx, K_m0, p_2U, K, alpha, beta, gamma, k_e1, k_e2
    V_G, k_1, k_2, V_I, m_1, m_2, m_4, m_5, m_6, HE_b, k_max, k_min, k_abs, k_gri, f, a, b, c, d, k_p1, k_p2, k_p3, k_p4, k_i, F_cns, V_m0, V_mx, K_m0, p_2U, K, alpha, beta, gamma, k_e1, k_e2 = return_parameters("norm")
    global G_pb, G_b, I_lb, I_pb, I_b, G_tb, I_pob
    G_pb = 4.5
    G_b = 4.5
    I_lb = 15
    I_pb = 15
    I_b = 15
    G_tb = 4.5
    I_pob = 15

    for i in range (13):
        print (G(i))
    
    V_G, k_1, k_2, V_I, m_1, m_2, m_4, m_5, m_6, HE_b, k_max, k_min, k_abs, k_gri, f, a, b, c, d, k_p1, k_p2, k_p3, k_p4, k_i, F_cns, V_m0, V_mx, K_m0, p_2U, K, alpha, beta, gamma, k_e1, k_e2 = return_parameters("diab")

    for i in range(13):
        print (G(i))
    import matplotlib.pyplot as plt

    iterations = range(5)
    glucose_levels = [G(i) for i in iterations]

    plt.plot(iterations, glucose_levels, marker='o')
    plt.xlabel('Iteration')
    plt.ylabel('Glucose Level')
    plt.title('Glucose Level Over Iterations')
    plt.grid(True)
    plt.show()
    plt.waitforbuttonpress()