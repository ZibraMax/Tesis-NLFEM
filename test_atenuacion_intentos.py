import numpy as np
import matplotlib.pyplot as plt

N = 50
A = 1

l = 0.1
t = 0.5
L0 = 1/2/np.pi/l**2

_x = np.linspace(-A,A,N)
_y = np.linspace(-A,A,N)

X = []
Y = []
Z = []
Z1 = []
Z2 = []

r = lambda x,y: np.sqrt(x**2+y**2)
macauley = lambda f, *args: f(*args) + abs(f(*args))

f = lambda x,y: L0*r(x,y)*np.exp(-r(x,y)/l)
f1 = lambda x,y: L0*np.exp(-r(x,y)/l)

f2p = lambda x,y: 1-(2*l*r(x,y)/l)
f2 = lambda x,y: macauley(f2p,x,y)

for x in _x:
	for y in _y:
		X.append(x)
		Y.append(y)
		Z.append(f(x,y))
		Z1.append(f1(x,y))
		Z2.append(f2(x,y))
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# surf = ax.plot_trisurf(X, Y, Z,cmap='magma')
# surf = ax.plot_trisurf(X, Y, Z1,cmap='magma')
surf = ax.plot_trisurf(X, Y, Z2,cmap='magma')
surf._facecolors2d=surf._facecolors3d
surf._edgecolors2d=surf._edgecolors3d
cbar = fig.colorbar(surf)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title(r'$f(x)=\frac{1}{2\pi l^2 t}\left(\sqrt{x^2+y^2}\right)e^{\left(\frac{-\sqrt{x^2+y^2}}{l}\right)}$')
plt.show()

# psis = lambda z,n: np.array([[1/4*(1-z)*(1-n)*(-1-z-n)],[1/4*(1+z)*(1-n)*(-1+z-n)],[1/4*(1+z)*(1+n)*(-1+z+n)],[1/4*(1-z)*(1+n)*(-1-z+n)],[1/2*(1-z**2)*(1-n)],[1/2*(1+z)*(1-n**2)],[1/2*(1-z**2)*(1+n)],[1/2*(1-z)*(1-n**2)]])
# dzpsis = lambda z,n: np.array([[-1/4*(n-1)*(2*z+n)],[-1/4*(n-1)*(2*z-n)],[1/4*(n+1)*(2*z+n)],[1/4*(n+1)*(2*z-n)],[(n-1)*z],[-1/2*(n**2-1)],[-(n+1)*z],[1/2*(n**2-1)]])
# dnpsis = lambda z,n: np.array([[-1/4*(z-1)*(2*n+z)],[1/4*(z+1)*(2*n-z)],[1/4*(z+1)*(2*n+z)],[-1/4*(z-1)*(2*n-z)],[1/2*(z**2-1)],[-n*(z+1)],[-1/2*(z**2-1)],[n*(z-1)]])
# LR = 6*l
# dl = np.sqrt(2)/2*LR
# coords = np.array([
# 	[-dl,-dl],
# 	[dl,-dl],
# 	[dl,dl],
# 	[-dl,dl],
# 	[0,-LR],
# 	[LR,0],
# 	[0,LR],
# 	[-LR,0]])

# dxdz = lambda z, n: coords[:, 0] @ dzpsis(z, n)
# dydz = lambda z, n: coords[:, 1] @ dzpsis(z, n)
# dxdn = lambda z, n: coords[:, 0] @ dnpsis(z, n)
# dydn = lambda z, n: coords[:, 1] @ dnpsis(z, n)
# J = lambda z, n: np.array([[dxdz(z, n)[0], dydz(z, n)[0]], [dxdn(z, n)[0], dydn(z, n)[0]]])
# TX = lambda z, n: coords[:, 0] @ psis(z, n)
# TY = lambda z, n: coords[:, 1] @ psis(z, n)

# integral = 0
# n = 80
# P,W = np.polynomial.legendre.leggauss(n)

# for i in range(n):
# 	for j in range(n):
# 		integral+=f(TX(P[i],P[j])[0],TY(P[i],P[j])[0])*W[i]*W[j]*np.linalg.det(J(P[i],P[j]))

# print(integral)
