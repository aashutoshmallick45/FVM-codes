import numpy as np
import matplotlib.pyplot as plt

#length of the domain
l=1

#mesh
n=25
dx=l/(n-1)

#constants
rho = 1
u = 0.1
gamma = 0.1
F=rho * u
D=gamma/dx
aw=D+F/2
ae=D-F/2

#phi vales or variable
phi=np.zeros(n)
phiold = np.zeros(n)

# boundary conditions
phi[0]= 1
phi[n-1] = 0

su=0
sp=0
errorhistory = []
error = 1
iterationhistory = []
iteration = 0

while error > 1e-6:
    phitemp = phi.copy()
    for i in range(1,n-1):
        if i==1 :
           # aw = 0
            su = (2*D+F)*phi[0]
            sp = -(2*D+F)
            ap = ae - sp
            phi[i] = (ae*phi[i+1]+ su)/ap    # aw*phi[i-1] = 0 
        
        elif i== (n-2):

           # ae = 0
            su = (2*D-F)*phi[n-1]
            sp = -(2*D-F)
            ap = aw - sp
            phi[i] = (aw*phi[i-1] + su)/ap    # ae*phi[i+1] = 0

        else:
            su = 0
            sp = 0
            ap = aw + ae - sp
            phi[i] = (aw*phi[i-1] + ae*phi[i+1]+ su)/ap
            

    error = np.max(np.abs(phi - phitemp))
    errorhistory.append(error)
    iteration = iteration + 1
    iterationhistory.append(iteration)
    
    phiold=phi.copy()
print("soluciton from FVM iteration : ",phi)

#Analytical result
x_analytical = np.linspace(0, l, n-2) # Use more points for a smooth curve

# The formula
phi_exact = 1 - (np.exp(rho * u * x_analytical / gamma) - 1) / (np.exp(rho * u * l / gamma) - 1)

print("Analytical solution is : ", phi_exact)

# --- Visualization ---

# Plot 1: Convergence (Residuals)
plt.figure(1)
plt.semilogy(iterationhistory, errorhistory, color='blue') # Log scale is best for residuals
plt.title('Convergence ')
plt.xlabel('Iterations')
plt.ylabel('Max Residual Error')
plt.grid(True)

# Plot 2: Final Result (Phi vs x)
plt.figure(2)
x_coords = np.linspace(0, l, n)
plt.plot(x_coords, phi, 'o-r', label='FVM Result')
plt.title('FVM Steady-State Solution for u= 0.1m/s')
plt.xlabel('Distance x')
plt.ylabel('Phi')
plt.legend()
plt.grid(True)

plt.figure(3)
plt.plot(x_coords, phi, 'o-r', label='From FVM')
plt.plot(x_analytical, phi_exact, 'g-', label='Analytical result')
plt.title('Comparision Analytical vs FVM result at u = 0.1m/s')
plt.xlabel('Distance x')
plt.ylabel('Phi')
plt.legend()
plt.grid(True)

plt.show()


