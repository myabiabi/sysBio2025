import random
import math

# ========================
# Parámetros
# ========================
k1 = 2.0   # tasa de producción (moleculas/min)
k2 = 0.1   # tasa de degradación (1/min)
NA = 0     # moléculas iniciales
t = 0.0    # tiempo inicial
t_max = 50 # tiempo máximo de simulación

# Guardar trayectoria
tiempos = [t]
moleculas = [NA]

# ========================
# Algoritmo de Gillespie
# ========================
while t < t_max:
    # 1. Calcular propensiones
    lambda1 = k1
    lambda2 = k2 * NA
    lambdas = [lambda1, lambda2]
    lambda_T = sum(lambdas)

    if lambda_T <= 0:
        break  # no hay más reacciones posibles
    
    # 2. Tiempo hasta la siguiente reacción
    r1 = random.random()
    tau = (1.0 / lambda_T) * math.log(1.0 / r1)
    t += tau

    # 3. Elegir qué reacción ocurre
    r2 = random.random() * lambda_T
    if r2 < lambda1:
        # Reacción 1: Producción
        NA += 1
    else:
        # Reacción 2: Degradación
        if NA > 0:
            NA -= 1

    # 4. Guardar datos
    tiempos.append(t)
    moleculas.append(NA)

# ========================
# Mostrar resultados
# ========================
for ti, ni in zip(tiempos, moleculas):
    print(f"t={ti:.3f}, A={ni}")



##############
import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random

# ========================
# Parámetros del sistema
# ========================
k = 2.0       # tasa de producción
gamma = 0.1   # tasa de degradación
X = [0]       # moléculas iniciales de X
t = [0.0]     # tiempo inicial
tend = 1000   # tiempo máximo de simulación

# ========================
# Algoritmo de Gillespie
# ========================
while t[-1] < tend:
    current_X = X[-1]
    
    # 1. Calcular propensiones
    rates = [k, gamma * current_X]
    rate_sum = sum(rates)
    
    if rate_sum <= 0:
        break  # no hay más reacciones posibles
    
    # 2. Tiempo hasta la siguiente reacción (exponencial con tasa=rate_sum)
    tau = np.random.exponential(scale=1/rate_sum)
    t.append(t[-1] + tau)
    
    # 3. Elegir qué reacción ocurre
    r2 = random.uniform(0, 1) * rate_sum
    
    if r2 <= rates[0]:
        # Reacción de producción
        X.append(current_X + 1)
    else:
        # Reacción de degradación
        X.append(current_X - 1)

# ========================
# Resultados
# ========================
plt.plot(t, X, drawstyle="steps-post")
plt.xlabel("Tiempo")
plt.ylabel("Cantidad de X (mRNA)")
plt.title("Simulación Gillespie (nacimiento-muerte)")
plt.show()
