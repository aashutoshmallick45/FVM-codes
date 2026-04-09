import numpy as np
import matplotlib.pyplot as plt

l = 0.5    
n = 50  
dx = l / (n - 1)

# Material properties
k = 1000
A = 10 * 10**-3

# Boundary conditions
Ta = 100
Tb = 500

T = np.zeros(n)
T[0] = Ta
T[-1] = Tb # Uses the last index automatically
error = 1

ae = k * A / dx
aw = k * A / dx
ap = aw + ae

# Solver
while error > 1e-6:
    Told = T.copy()
    for i in range(1, n - 1): # Dynamic range
        T[i] = (aw * T[i-1] + ae * T[i+1]) / ap
    error = np.max(np.abs(T - Told))

# Plotting - Fixed X axis
X = np.linspace(0, l, n) # 5 points from 0 to 0.5
plt.plot(X, T, marker='o', color='red', label='Numerical Solution')
plt.xlabel('Position along rod (m)')
plt.ylabel('Temperature (°C)')
plt.title('1D Steady State Heat Conduction')
plt.grid(True)
plt.legend()
plt.show()

print("Final Temperatures:", T)



