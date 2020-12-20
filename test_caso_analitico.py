import numpy as np
import matplotlib.pyplot as plt

L = 50
E = 2.1*10**6
l = 0.1
Z1 = 0.5
sbarra = 2100

def ex_analitico(x,L,E,l,Z1,sbarra):
	l0 = 1/2/l
	ebarra = sbarra/E
	ld = -l0*(1-Z1)/Z1
	return ebarra-ld*l/2*ebarra*(np.exp((ld*x*l-x)/l)+np.exp((ld*l*L-ld*l*x-L+x)/l))
_X = np.linspace(0,L,1000)
plt.plot(_X,ex_analitico(_X,L,E,l,Z1,sbarra))
plt.grid()
# plt.ylim([4.5*10**-5,6.5*10**-5])
# plt.savefig('defx_analitico.png',transparent=True)
plt.show()
