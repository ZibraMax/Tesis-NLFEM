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
f1 = lambda rho: np.exp(-rho)
f2p = lambda rho: (1-(rho/(m0)))
f3p = lambda rho: (1-(rho**2/(m0)**2))
f2 = lambda rho: macauley(f2p,rho)
f3 = lambda rho: macauley(f3p,rho)
f4 = lambda rho: rho*np.exp(-rho)

L0S = np.array([L0,L0_lineal,L0_cuadratico,L0_custom])
atenuacion = lambda rho: np.array([f1(rho),f2(rho),f3(rho),f4(rho)])*L0S
# ls = [0.1]
ls = [0.1,0.2,0.3,0.4,0.5]
fig = plt.figure()
ax = fig.add_subplot()
Z = np.linspace(-8,8,1000)
for i in range(len(atenuacion(0))):
	for l in ls:
		m0 = 6
		LR = (m0)*l
		L0 = 1/2/np.pi/l/l
		L0_lineal = 3/2/np.pi/LR/LR
		L0_cuadratico = 1/np.pi/LR/LR
		L0_custom = 1/4/np.pi/l/l

		r = lambda x,y: np.sqrt(x**2+y**2)
		macauley = lambda f, *args: f(*args)/2 + abs(f(*args))/2
		f1 = lambda rho: np.exp(-rho)
		f2p = lambda rho: (1-(rho/(m0)))
		f3p = lambda rho: (1-(rho**2/(m0)**2))
		f2 = lambda rho: macauley(f2p,rho)
		f3 = lambda rho: macauley(f3p,rho)
		f4 = lambda rho: rho*np.exp(-rho)

		# L0S = np.array([1,1,1,1])
		L0S = np.array([L0,L0_lineal,L0_cuadratico,L0_custom])
		atenuacion = lambda rho: np.array([f1(rho),f2(rho),f3(rho),f4(rho)])*L0S
		X = []
		Y = []
		for z in Z:
			X.append(z)
			Y.append(atenuacion(abs(z))[i].tolist())
		plt.plot(X,Y,label=r'$l=$'+format(l))
	plt.legend()
	plt.xlabel(r"$\rho$")
	plt.ylabel(r"$A(|x-x'|)$")
	plt.grid()
	plt.show()