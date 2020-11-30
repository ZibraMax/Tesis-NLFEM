import numpy as np
import matplotlib.pyplot as plt


m0 = 6
l = 0.5
LR = (m0)*l
L0 = 1/2/np.pi/l/l
L0_lineal = 3/2/np.pi/LR/LR
L0_cuadratico = 1/np.pi/LR/LR
L0_custom = 1/4/np.pi/l/l


r = lambda x,y: np.sqrt(x**2+y**2)
macauley = lambda f, *args: f(*args) + abs(f(*args))
f1 = lambda x,y: np.exp(-r(x,y)/l)
f2p = lambda x,y: (1-(r(x,y)/LR))
f3p = lambda x,y: (1-(r(x,y)**2/LR**2))
f2 = lambda x,y: macauley(f2p,x,y)
f3 = lambda x,y: macauley(f3p,x,y)
f4 = lambda x,y: r(x,y)/l*np.exp(-r(x,y)/l)

L0S = np.array([L0,L0_lineal,L0_cuadratico,L0_custom])
atenuacion = lambda x,y: np.array([f1(x,y),f2(x,y),f3(x,y),f4(x,y)])*L0S

psis = lambda z,n: np.array([[1/4*(1-z)*(1-n)*(-1-z-n)],[1/4*(1+z)*(1-n)*(-1+z-n)],[1/4*(1+z)*(1+n)*(-1+z+n)],[1/4*(1-z)*(1+n)*(-1-z+n)],[1/2*(1-z**2)*(1-n)],[1/2*(1+z)*(1-n**2)],[1/2*(1-z**2)*(1+n)],[1/2*(1-z)*(1-n**2)]])
dzpsis = lambda z,n: np.array([[-1/4*(n-1)*(2*z+n)],[-1/4*(n-1)*(2*z-n)],[1/4*(n+1)*(2*z+n)],[1/4*(n+1)*(2*z-n)],[(n-1)*z],[-1/2*(n**2-1)],[-(n+1)*z],[1/2*(n**2-1)]])
dnpsis = lambda z,n: np.array([[-1/4*(z-1)*(2*n+z)],[1/4*(z+1)*(2*n-z)],[1/4*(z+1)*(2*n+z)],[-1/4*(z-1)*(2*n-z)],[1/2*(z**2-1)],[-n*(z+1)],[-1/2*(z**2-1)],[n*(z-1)]])
dl = np.sqrt(2)/2*LR
coords = np.array([
	[-dl,-dl],
	[dl,-dl],
	[dl,dl],
	[-dl,dl],
	[0,-LR],
	[LR,0],
	[0,LR],
	[-LR,0]])

dxdz = lambda z, n: coords[:, 0] @ dzpsis(z, n)
dydz = lambda z, n: coords[:, 1] @ dzpsis(z, n)
dxdn = lambda z, n: coords[:, 0] @ dnpsis(z, n)
dydn = lambda z, n: coords[:, 1] @ dnpsis(z, n)
J = lambda z, n: np.array([[dxdz(z, n)[0], dydz(z, n)[0]], [dxdn(z, n)[0], dydn(z, n)[0]]])
TX = lambda z, n: coords[:, 0] @ psis(z, n)
TY = lambda z, n: coords[:, 1] @ psis(z, n)
n_atenuacion = len(atenuacion(0,0))
integral = np.zeros([n_atenuacion])
n = 10*m0
print('Integrando con ' + format(n) + ' puntos de Gauss')
P,W = np.polynomial.legendre.leggauss(n)
for i in range(n):
	for j in range(n):
		integral += atenuacion(TX(P[i],P[j])[0],TY(P[i],P[j])[0])*W[i]*W[j]*np.linalg.det(J(P[i],P[j]))
print('INTEGRALES:')
print(integral)
print('')
print("L0's:")
print(L0S)

n_puntos = 300
Z = np.linspace(-1,1,n_puntos)
N = np.linspace(-1,1,n_puntos)

X = []
Y = []
U = []

for z in Z:
	for n in N:
		X.append(TX(z,n)[0])
		Y.append(TY(z,n)[0])
		U.append(atenuacion(X[-1],Y[-1]).tolist())
U = np.array(U)
print("min(f(x,y))")
print(np.min(U,axis=0))

lado = int(np.ceil(np.sqrt(n_atenuacion)))
fig = plt.figure()
for i in range(n_atenuacion):
	ax = fig.add_subplot(lado,lado,i+1,projection='3d')
	surf = ax.plot_trisurf(X,Y,U[:,i],cmap='magma')
	ax.set_title(r"$\int_{V'}{A(|x-x'|)}{dv'}=$"+'{:.3f}'.format(integral[i])+'\n'+r'$\lambda_0=$'+'{:.3f}'.format(L0S[i]))
	plt.colorbar(surf)
	ax.set_xlabel(r'$x$')
	ax.set_ylabel(r'$y$')
	ax.set_zlabel(r'$f(x,y)$')
fig.suptitle('Funciones de Atenuación con:\n'+r'$l=$'+format(l)+'\n'+r'$LR='+format(m0)+'l$')
plt.show()
