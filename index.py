import NLFEM
import numpy as np
import matplotlib.pyplot as plt
from NLFEM.Mesh import delaunay


def generarCBdesdeBordeX(this, borde, valor=0):
    cb = []
    nodos = this.darNodosCB(borde)
    cbe = np.zeros([len(nodos), 2])
    cbe[:, 0] = nodos*2
    cbe[:, 1] = valor
    cb += cbe.tolist()
    return cb

def generarCBdesdeBordeY(this, borde, valor=0):
    cb = []
    nodos = this.darNodosCB(borde)
    cbe = np.zeros([len(nodos), 2])
    cbe[:, 0] = nodos*2+1
    cbe[:, 1] = valor
    cb += cbe.tolist()
    return cb

def generarCBdesdeBorde(this, borde, valor=[0,0]):
    return generarCBdesdeBordeX(this, borde, valor[0])+generarCBdesdeBordeY(this, borde, valor[1])


a = 5
E = 21*10**6
V = 0.2
u = 0.001

GEOMETRIA = delaunay.Delaunay1V([[0,0],[a,0],[a,a],[0,a]],delaunay._strdelaunay(a=4,o=2),plot=True)
plt.show()
Objeto_FEM = NLFEM.NoLocal(GEOMETRIA)
Objeto_FEM.generarElementos()
GEOMETRIA.cbe = generarCBdesdeBorde(GEOMETRIA,3,[0,0])+generarCBdesdeBorde(GEOMETRIA,1,[u,0])
Objeto_FEM.definirCondicionesDeBorde(GEOMETRIA.cbe)
Objeto_FEM.solucionar(E=E,v=V,Fx=lambda x,y: 0,Fy=lambda x,y: 0,plot=False)
z1=0.6
Ke = Objeto_FEM.elementos[0].Ke*z1+Objeto_FEM.elementos[0].KNLS[0]*(1-z1)

#Objeto_FEM.defUnitariaX()
plt.show()
print(Ke)
a=a