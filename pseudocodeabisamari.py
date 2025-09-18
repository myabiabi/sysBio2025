import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random
import math
# =========================
# Parámetros
# donde k1 es la tasa en la que se produce mRNA y k2 es la tasa en la que se degrada
# ==================
k1 = 2
k2 = 0.1
# ===================
# Moléculas y tiempo inicial
# tanto X_0 como t_0 se actualizarán con cada evento que pase
# ==================
mrna=0
t=0
MRNA = [0]
T = [0]
current_X = MRNA[-1]
# ==============
# tiempo de la simulación
# ==============
t_final = 1000
# ==============
# propensiones 
# ==============
a1 = k1
a2 = k2 * current_X
#Algoritmo de Gillespie
# =====================
while T[-1] < t_final: # el indice -1 llama al utlimo elemento del vector t, el cual al inicio es 0
        current_X = MRNA[-1] #current_X se refiere al número de moléculas del útlimo evento, también la llamamos con el indice -1
#1. Calcular propensiones
#propension de la reaccion de producción de mRNA y de degradación de mRNA
#propension total, suma de las dos anteriores
#si la propension total es menor o igual a 0, ninguna reacción puede ocurrir más, si k es 0 no se produce más mRNA,
#si X es cero, es porque ya no hay moleculas que reaccionen 
        a_sigma = a1 + a2  
        if(a_sigma <= 0):
                break 
#2. tiempo de reacción 
#el tiempo de espera de los entre evento y evento sigue una distribución esponencial

        r1 = random.uniform(0,1)
        tau = -math.log(r1)/a_sigma
        T.append(T[-1] + tau)

#3. elegir cuál reacción ocurre

        r2 = random.uniform(0,1) * a_sigma
        if r2 <= a1:
            MRNA.append(current_X + 1)
        else:
            MRNA.append(current_X - 1)



plt.plot(T, MRNA, drawstyle="steps-post")
plt.xlabel("Tiempo")
plt.ylabel("Cantidad de X (mRNA)")
plt.title("Simulación Algoritmo de Gillespie")
plt.show()
