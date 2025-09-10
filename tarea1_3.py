import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ========================
# Sistema ODE con estímulo variable
# ========================
def sistema(t, x, kA_on, kB, delta, mu, KD, n, t_off, t_on):
    A, B, C = x

    # Estímulo externo (apagado entre t_off y t_on)
    if t_off <= t <= t_on:
        kA = 0
    else:
        kA = kA_on

    # Dinámica de A
    dAdt = kA - delta * A

    # Dinámica de B
    dBdt = kB * A - delta * B

    # Dinámica de C (Hill con A y/o B)
    act = (A**n + B**n) / (A**n + B**n + KD**n)
    dCdt = mu * act - delta * C

    return [dAdt, dBdt, dCdt]


# ========================
# Parámetros
# ========================
delta = 0.1          # tasa de dilución (1/min)
kA_on = delta * 20   # activa A → 20 nM
kB = 0.1
KD = 10
n = 16
mu = delta * 15      # saturación en 15 nM

# Condiciones iniciales
x0 = [0, 0, 0]

# Tiempo de simulación
t_eval = np.linspace(0, 400, 1000)

# ========================
# Simulación 1: estímulo OFF por 10 min
# ========================
sol1 = solve_ivp(sistema, [0, max(t_eval)], x0,
                 args=(kA_on, kB, delta, mu, KD, n, 100, 110),  # OFF entre 100 y 110 min
                 t_eval=t_eval)

# ========================
# Simulación 2: estímulo OFF por 60 min
# ========================
sol2 = solve_ivp(sistema, [0, max(t_eval)], x0,
                 args=(kA_on, kB, delta, mu, KD, n, 100, 160),  # OFF entre 100 y 160 min
                 t_eval=t_eval)

# ========================
# Gráficas
# ========================
fig, axs = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

for sol, label in zip([sol1, sol2], ["OFF 10 min", "OFF 60 min"]):
    t, (A, B, C) = sol.t, sol.y
    axs[0].plot(t, A, lw=2, label=f"A ({label})")
    axs[0].plot(t, B, lw=2, label=f"B ({label})")
    axs[0].plot(t, C, lw=2, label=f"C ({label})")

axs[0].set_ylabel("Concentración (nM)")
axs[0].legend()
axs[0].grid(True, linestyle="--", alpha=0.6)
axs[0].set_title("Respuesta de A, B y C con apagado del estímulo")

# Señal externa (para visualizar cuándo se apaga)
stimulus = np.ones_like(t_eval) * kA_on
stimulus[(t_eval >= 100) & (t_eval <= 110)] = 0  # caso 1
axs[1].plot(t_eval, stimulus, label="Estímulo (OFF 10 min)", lw=2)

stimulus2 = np.ones_like(t_eval) * kA_on
stimulus2[(t_eval >= 100) & (t_eval <= 160)] = 0  # caso 2
axs[1].plot(t_eval, stimulus2, label="Estímulo (OFF 60 min)", lw=2)

axs[1].set_xlabel("Tiempo (min)")
axs[1].set_ylabel("kA (nM/min)")
axs[1].legend()
axs[1].grid(True, linestyle="--", alpha=0.6)

plt.tight_layout()
plt.show()

#