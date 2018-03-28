# Created by Matthew Herald, mherald1@vols.utk.edu
# Solves one-group neutron transport equation for non-multiplying medium in a 1-D sphere
# User input is required for the number of nodes and sphere radius
# Analytical and numerical solutions are generated for comparison
# 3-28-18

import numpy as np
import matplotlib.pyplot as plt

nodes = int(input('Input number of nodes: '))
R = float(input('Input radius of sphere (cm): '))

macro_ab = 0.1532                   # cross-sections
macro_tr = 0.0362
macro_tot = macro_ab+macro_tr
D = 1/(3*macro_tot)                 # diffusion coefficient
L = np.sqrt(D/macro_ab)             # diffusion length
dr = R/(nodes-1)                    # delta r
r = np.zeros([1, nodes])

for i in range(0, nodes):           # creates array of node positions
    r[0, i] = dr*i

s = 1*10**10                        # point source at center of sphere
sourceTerm = np.zeros([nodes-1, 1]) # source terms array
sourceTerm[0] = s/(4*np.pi)         # source term at the center
mat = np.zeros([(nodes-1), (nodes-1)])  # creates matrix for the system of equations

for i in range(0, (nodes-1)):       # down

    n1 = (-D * ((r[0, i] - (dr / 2)) ** 2) / dr)  # interior node term 1
    n2 = ((D * ((r[0, i] + (dr / 2)) ** 2) / dr) + (D * ((r[0, i] - (dr / 2)) ** 2) / dr)
          + (macro_ab / 3) * ((r[0, i] + (dr / 2)) ** 3 - (r[0, i] - (dr / 2)) ** 3))
    n3 = (-D * ((r[0, i] + (dr / 2)) ** 2) / dr)  # interior node term 3

    for j in range(0, nodes-1):     # across
        if i == 0:                  # central node matrix values
            if j == 0:
                mat[i, j] = (((macro_ab*(dr**3))/24)+(D*dr/4))  # central node input
            elif j == 1:
                mat[i, j] = (-D*dr/4)    # central node input 2
            else:
                mat[i, j] = 0       # all other matrix values = 0

        else:
            if j == (i-1):
                mat[i, j] = n1      # interior node matrix value 1
            elif j == i:
                mat[i, j] = n2      # interior node matrix value 2
            elif j == (i+1):
                mat[i, j] = n3      # interior node matrix value 3
            else:
                mat[i, j] = 0       # all other matrix values = 0

inv = np.linalg.inv(mat)            # inverts matrix
numeric = inv*sourceTerm            # solves for the flux at each node

numeric = numeric[0, :]             # put matrix in format for plotting


def analytical(radius):

    # solves for the analytical solution
    a = (-s*np.exp(2*R/L)/(4*np.pi*D*(1-np.exp(2*R/L))))
    b = s/(4*np.pi*D*(1-np.exp(2*R)))
    flux = a*np.exp(-radius/L)/r + b*np.exp(radius/L)/radius

    return flux


divs = 100  # number of points to be evaluated
x = np.linspace((1/50)*R, R, divs)  # analytical sphere radius points
y = analytical(x)                   # analytical flux


# numerical solution plot
p1 = plt.scatter(r[0, 0:len(r[0]) - 1], numeric, marker='*', c='r')
# analytical solution plot
p2 = plt.scatter(x, y, marker='o', c='b')
plt.title('Numeric and Analytical Solutions to Neutron Diffusion')
plt.ylabel('Neutron flux (#n/s)')
plt.xlabel('Node Position (cm)')
plt.legend([p1, p2], ['Numerical', 'Analytical'])
plt.yscale('log')   # plots the y-axis on a log scale
plt.show()