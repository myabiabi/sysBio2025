import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the Lotka-Volterra equations
def lotka_volterra(X, t, alpha1, theta1, alpha2, theta2, alpha3, theta3):
    B, V = X  # Unpack prey and predator populations

    # Define the differential equations
    dBdt = alpha1 * B - theta1 * B * V  # Prey population change
    dVdt = alpha2 * V + theta2 * B * V  # Predator population change

    return [dBdt, dVdt]

# Set model parameters (example values)
alpha1 = 0.15  # Prey birth rate
alpha2 = 0.15
alpha3 = 0.15
theta1 = 0.30
theta2 = 0.13
theta3 = 0.03

# Set initial conditions
B0 = 10.0  # Initial prey population
V0 = 2.0   # Initial predator population
initial_conditions = [B0, V0]

# Define time points for simulation
t = np.linspace(0, 200, 1000)  # Simulate for 200 time units, with 1000 steps

# Solve the differential equations using odeint
# args=(alpha, beta, gamma, delta) passes these parameters to the lotka_volterra function
solution = odeint(lotka_volterra, initial_conditions, t, args=(alpha1,alpha2, alpha3, theta1, theta2, theta3))

# Extract prey and predator populations from the solution
prey_population = solution[:, 0]
predator_population = solution[:, 1]

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(t, prey_population, label='Bacteriofagos (B)')
plt.plot(t, predator_population, label='Virus (V)')
plt.xlabel('Tiempo')
plt.ylabel('Poblacion')
plt.title('Lotka-Volterra Bacteriofagos')
plt.legend()
plt.grid(True)
plt.show()
