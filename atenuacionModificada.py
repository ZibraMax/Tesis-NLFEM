import numpy as np
import matplotlib.pyplot as plt


m0 = 6
l = 0.1
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

L0S = np.array([L0_custom])
atenuacion = lambda x,y: np.array([f4(x,y)])*L0S

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

ls = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
ls2 = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
fig = plt.figure()
ax = fig.add_subplot()
Z = np.linspace(-1,1,1000)
for l,l2 in zip(ls,ls2):
	f4 = lambda x,y: r(x,y)/l*np.exp(-r(x,y)/l)
	L0_custom = 1/4/np.pi/l/l
	L0S = np.array([L0_custom])
	atenuacion = lambda x,y: np.array([f4(x,y)])*L0S
	LR = (m0)*l
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
	X = []
	Y = []
	for z in Z:
		X.append(TX(z,0)[0])
		Y.append(atenuacion(X[-1],0).tolist())
	plt.plot(X,Y,label=format(l)+'')
plt.legend()
plt.show()