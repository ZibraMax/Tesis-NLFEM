import NLFEM
import numpy as np
import matplotlib.pyplot as plt
x0=0
y0=0
a=2
b=2


i=1
j=1

coords = [[x0,y0],[x0+a,y0],[x0+a,y0+b],[x0,y0+b],[x0+a/2,y0],[x0+a,y0+b/2],[x0+a/2,y0+b],[x0,y0+b/2]]
elemento = NLFEM.SerendipityC(coords)
fig = plt.figure(figsize=[30,30])
f = lambda z,n,e: e._J(z,n)[0][0]*e.dzpsis(z,n)[i][0]+e._J(z,n)[0][1]*e.dnpsis(z,n)[i][0]
Z = elemento._dominioNaturalZ
N = elemento._dominioNaturalN
U = []
for z, n in zip(Z, N):
	U.append(f(z,n,elemento))
ax = fig.add_subplot(projection='3d')
surf = ax.plot_trisurf(Z, N, U, cmap='magma')
cbar = fig.colorbar(surf)
ax.set_xlabel(r'$\zeta$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$u$')
plt.show()