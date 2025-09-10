import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ========================
# Sistema ODE Caso #4
# ========================
def sistema(t, x, kA, kB, delta, mu, KD, n):
    A, B, C = x

    # Dinámica de A (constitutiva con estímulo externo)
    dAdt = kA - delta * A

    # Dinámica de B (lineal con A)
    dBdt = kB * A - delta * B

    # Dinámica de C (Hill cooperativo con A y/o B)
    act = (A**n + B**n) / (A**n + B**n + KD**n)
    dCdt = mu * act - delta * C

    return [dAdt, dBdt, dCdt]


# ========================
# Parámetros
# ========================
delta = 0.1        # tasa de dilución (1/min)
kA = delta * 20    # para que A -> 20 nM en equilibrio
kB = 0.1           # tasa de producción de B (1/min)
KD = 10            # nM
n = 16             # coeficiente de Hill
mu = delta * 15    # para que C -> 15 nM máximo

# Condiciones iniciales
A0 = 0
B0 = 0
C0 = 0
x0 = [A0, B0, C0]

# Tiempo de simulación
t_eval = np.linspace(0, 200, 500)

# Resolver ODE
sol = solve_ivp(sistema, [0, max(t_eval)], x0,
                args=(kA, kB, delta, mu, KD, n),
                t_eval=t_eval)

# ========================
# Resultados
# ========================
t = sol.t
A, B, C = sol.y

# Gráficas
plt.figure(figsize=(8,5))
plt.plot(t, A, label="A", lw=2)
plt.plot(t, B, label="B", lw=2)
plt.plot(t, C, label="C", lw=2)
plt.xlabel("Tiempo (min)")
plt.ylabel("Concentración (nM)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.title("Caso #4: Dinámica de A, B y C")
plt.show()


#este funciona y grafica
