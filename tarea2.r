k1 <- 2.0 #tasa de produccion
k2 <- 0.1 #tasa de degradación
NO <- 0 #moleculas iniciales
t <- 0.0 #tiempo inicial
t_max <- 50 #tiempo maximo de simulacion

tiempo <- c(t)
moleculas <- c(NO)

while (t < t_max) {
    lamb1 <- k1
    lamb2 <- k2 * NO
    lamb_t <- lamb1 + lamb2 

    if (lamb_t <= 0){
        break
    }

}

r1 <- runif(1)
tau <- (1 / lamb_t) * log(1 / r1)
t <- t + tau

r2 <- runif(1) * lambda_t
if (r2 < lamb1){
    N0 <- NO + 1 
} else {
   if (NO > 0) {
      NO <- NO -1
   }
}

tiempos <- c (tiempo, t)
moleculas <- c(moleculas, NO)

resultado <- data.frame(tiempo = tiempos, A = moleculas)
print(resultado)
