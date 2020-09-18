import NLFEM
import numpy as np
import matplotlib.pyplot as plt
from NLFEM.Mesh import Rect
import os

PATH = os.getcwd()
nombre = "\\NLFEM\\Mesh\\input.txt"
GEOMETRIA = Rect.Rect(PATH+nombre)

GEOMETRIA.dibujarse()

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
E = 2.1*10**6
V = 0.2
u = 0.001

Objeto_FEM = NLFEM.NoLocal(GEOMETRIA)
GEOMETRIA.cbe = generarCBdesdeBorde(GEOMETRIA,3,[0,0])+generarCBdesdeBordeX(GEOMETRIA,1,u)
Objeto_FEM.definirCondicionesDeBorde(GEOMETRIA.cbe)
z1=0.6
Objeto_FEM.z1 = z1
Objeto_FEM.solucionar(E=E,v=V,Fx=lambda x,y: 0,Fy=lambda x,y: 0,plot=True)
Ke = Objeto_FEM.elementos[0].Ke*z1+Objeto_FEM.elementos[0].KNLS[0]*(1-z1)
Objeto_FEM.defUnitariaX()
plt.show()