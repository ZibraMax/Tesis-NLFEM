from NLFEM.FEM1V import *
import numpy as np
import matplotlib.pyplot as plt
from NLFEM.Mesh import Rect
import os

PATH = os.getcwd()
nombre = "\\NLFEM\\Mesh\\input.txt"
GEOMETRIA = Rect.Rect(PATH+nombre)

GEOMETRIA.dibujarse()

GEOMETRIA.definirTodasCondiciones()
zanahorias = FEM1V(GEOMETRIA)
zanahorias.definirCondicionesDeBorde(GEOMETRIA.cbe)
E = 200000
v = 0.27
G = E/(2*(1+v))
a11 = lambda x,y: 1
a12 = lambda x,y: 0
a21 = lambda x,y: 0
a22 = lambda x,y: 1
a00 = lambda x,y: 0

theta = 1
f = lambda x,y: 2*G*theta

zanahorias.solucionar(cmap='magma',markersize=1,linewidth=1,a11=a11,a12=a12,a21=a21,a22=a22,a00=a00,f=f,plot=True)
plt.show()
