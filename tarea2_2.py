import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import random
import math
# =========================
# Parámetros
# donde k es la tasa en la que se produce mRNA y gamma es la tasa en la que se degrada
# ==================
k = 2
gamma = 0.1
# ===================
# Moleculas y tiempo inicial
# tanto X_0 como t_0 se actualizarán con cada evento que pase
# ==================
X_0 = [0]
t_0 = [0]
current_X = X_0[-1]
# ==============
# tiempo de la simulación
# ==============
t_max = 1000
# ==============
# propensiones 
# ==============
lambda1 = k
lambda2 = gamma * current_X
#Algoritmo de Gillespe
# =====================
while t_0[-1] < t_max: #el indice -1 llama al utlimo elemento del vector t, el cual al inicio es 0
        current_X = X_0[-1] #current_X se refiere al número de moléculas del útlimo evento, también la llamamos con el indice -1
#1. Calcular propensiones
#propension de la reaccion de producción de mRNA y de degradación de mRNA
#propension total, suma de las dos anteriores
#si la propension total es menor o igual a 0, ninguna reacción puede ocurrir más, si k es 0 no se produce más mRNA,
#si X es cero, es porque ya no hay moleculas que reaccionen 
        rates = [lambda1, gamma * current_X] 
        rate_sum = sum(rates)  
        if(rate_sum <= 0):
                break 
#2. tiempo de reacción 
#el tiempo de espera de los entre evento y evento sigue una distribución esponencial
#el objeto r1 genera un número random entre 0 y 1 incluyendolos, lo cual representa un tiempo de reacción aleatorio
#este tiempo de reacción aleatorio sigue una distribución exponencial
# r1 = 1 - e^-λ_t*τ
# al ya conocer r1, que despejar τ que es el tiempo de reacción
#despeje:
# e^-λ_t*τ = 1 - r1
#aplicamos logoritmo porque es la operación inversa de la exponencial, para lograr despejar τ
# -λ_t*τ = ln(1-r1)
#depejamos τ
# τ = -1/λ_t*ln(1 - r1) 
#como 1-r1 es uniforme la puedo remplazar con r1'
# τ = -1/λ_t*ln(r1') 
#por propiedad de logaritmos
#-ln(x) = ln(1/x), podemos simplificar la expresión 
# τ = 1/λ_t*ln(1/r1') 
        r1 = random.random()
        tau = (1.0 / rate_sum) * math.log(1.0 / r1)
        t_0.append(t_0[-1] + tau)

#3. elegir que reacción ocurre
#random.uniform regresa un numero random entre dos numeros especificados, incluyendolos, esto se multiplica por la suma de la propensidades
#por lo que al final da un número de entre 0 y la sema de las propensidades
#rate[0] es para llamar el elemento 1 del obejto rate, en este caso lambda1 que es la propensidad de k
#si r2 es menor a igual que k, cae en el intervalo 0,k que representa la producción mientras que si cae en el intervalo de k, rate_sum que es degradación
        r2 = random.uniform(0,1) * rate_sum
        if r2 <= rates[0]:
            X_0.append(current_X + 1)
        else:
            X_0.append(current_X - 1)

plt.plot(t_0, X_0, drawstyle="steps-post")
plt.xlabel("Tiempo")
plt.ylabel("Cantidad de X (mRNA)")
plt.title("Simulación Algoritmo de Gillespie")
plt.show()
