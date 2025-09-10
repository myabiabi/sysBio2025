import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ========================
# Funciones auxiliares
# ========================
def fMM(x, mu, KD):
    """Michaelis-Menten like equation"""
    return mu * (x / (x + KD))

def fH(x, mu, KD, nH):
    """Hill like equation"""
    return mu * (x**nH / (x**nH + KD**nH)) #solo es conveniente definirla si la uso muchas veces 

def myODE_MM(t, x, k1, k2, k3, k4, n, mu, alpha):
    """Sistema de ecuaciones diferenciales para Michaelis-Menten"""
    A, B, C = x
    dAdt = (k1*E) - (k4*A)
    dBdt = (k2*A) - (k4*B)
    dCdt = mu * ( alpha + (1 - alpha) * ((A^n + B^n) / (A^n + B^n + KD^n)) ) - (k4*C) #aqui también esta la de hill
    return [dAdt, dBdt, dCdt]




# ========================
# Simulación ODE Michaelis-Menten
# ========================
k1 = 
k2 = 
k3 = 
k4 = 
n = 
mu= 
alpha =

S0 = 10000 # nM
E0 = 20    # nM
C0 = 0     # nM
x0 = [S0, E0, C0]

t_eval = np.arange(0, 10.1, 0.1)
sol = solve_ivp(myODE_MM, [0, 10], x0, args=(k1, k2, k3, k4, n, mu, alpha), t_eval=t_eval)

t = sol.t
A, B, C = sol.y
fMM_approx = mP * (E0 + C0) * (S0 + C0) / ((S0 + C0) + (kU / kB))

# Gráficas ODE vs Michaelis-Menten
fig, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

axs[0].plot(t, mP*C, '-', linewidth=2, color=(0.6, 0.8, 0), label="μP·C")
axs[0].plot(t, np.ones_like(t)*fMM_approx, ':', linewidth=2, color=(0.8, 0, 0.6), label="MM approx.")
axs[0].legend()
axs[0].set_ylabel("P synthesis")
axs[0].set_xlim([0, max(t)])
axs[0].grid(True, linestyle="--", alpha=0.5)

axs[1].plot(t, (fMM_approx - (mP*C)) / (mP*C), '--', linewidth=2, color=(0, 0.6, 0.8))
axs[1].set_xlabel("Time")
axs[1].set_ylabel("Relative error")
axs[1].set_xlim([0, max(t)])
axs[1].grid(True, linestyle="--", alpha=0.5)

plt.tight_layout()
plt.show()


# ========================
# Michaelis-Menten like function
# ========================
R = np.arange(0, 100.1, 0.1)
mu = 1
KD = 10

plt.figure(figsize=(6,4))
plt.plot(R, fMM(R, mu, KD), '-', linewidth=2, color=(0.6,0.8,0))
plt.plot(R, np.ones_like(R)*mu, '--', linewidth=2, color=(0,0.6,0.8))
plt.plot([0,KD,KD],[mu/2,mu/2,0],':', linewidth=2, color=(0.8,0,0.6))
plt.xlabel("Regulator [nM]")
plt.ylabel("f_MM")
plt.xlim([min(R), max(R)])
plt.ylim([0, mu*1.1])
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()


# ========================
# Hill function
# ========================
mu = 1
KD = 25
nH = 2

plt.figure(figsize=(6,4))
plt.plot(R, fH(R, mu, KD, nH), '-', linewidth=2, color=(0.6,0.8,0))
plt.plot(R, np.ones_like(R)*mu, '--', linewidth=2, color=(0,0.6,0.8))
plt.plot([0,KD,KD],[mu/2,mu/2,0],':', linewidth=2, color=(0.8,0,0.6))
plt.xlabel("Regulator [nM]")
plt.ylabel("f_H")
plt.xlim([min(R), max(R)])
plt.ylim([0, mu*1.1])
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()


# ========================
# Hill & Ultrasensitivity
# ========================
mu = 1
KD = 25
nH_vals = range(1,16,2)

plt.figure(figsize=(6,4))
colors = plt.cm.cool(np.linspace(0,1,len(nH_vals)))
for i, nH in enumerate(nH_vals):
    plt.plot(R, fH(R, mu, KD, nH), '-', linewidth=2, color=colors[i], label=f"nH={nH}")

plt.plot(R, np.ones_like(R)*mu, '--', linewidth=1, color=(0,0.6,0.8), label="μ")
plt.plot([0,KD,KD],[mu/2,mu/2,0],':', linewidth=1, color=(0.8,0,0.6), label="K_D")
plt.xlabel("Regulator [nM]")
plt.ylabel("f_H")
plt.xlim([min(R), max(R)])
plt.ylim([0, mu*1.1])
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()
