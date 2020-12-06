import NLFEM
import numpy as np
import matplotlib.pyplot as plt
from NLFEM.Mesh import Rect
import subprocess

RUTA_M = 'MATRICES'

E = 2.1*10**6
V = 0.2
u = 0.05
t = 0.1
l = 0.1

import os
if not os.path.exists(RUTA_M):
    subprocess.run("index.exe input.txt 3 "+format(E)+" "+format(V)+" "+format(t)+" "+format(l)+" "+RUTA_M+" 1", shell=True, check=True)

GEOMETRIA = Rect.Rect("input.txt")
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,0],[50,0])
GEOMETRIA.generarSegmentosDesdeCoordenadas([50,0],[50,1])
GEOMETRIA.generarSegmentosDesdeCoordenadas([50,1],[0,50])
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,50],[0,0])
# GEOMETRIA.dibujarse()
Objeto_FEM = NLFEM.NoLocal(GEOMETRIA)
condiciones_borde_escenciales = GEOMETRIA.generarCBdesdeBordeX(3,0)+GEOMETRIA.generarCBdesdeBordeX(1,u)+GEOMETRIA.generarCBYdesdeCoordenada(0,1/2)
Objeto_FEM.z1 = 0.5
Objeto_FEM.generarElementos()
Objeto_FEM.importarMatricesCarpeta(RUTA_M)
Objeto_FEM.definirCondicionesDeBorde(condiciones_borde_escenciales)
Objeto_FEM.ensamblar()
Objeto_FEM.condicionesFrontera()
Objeto_FEM.solucionarSistemaEcuaciones()

Objeto_FEM.defUnitariaX([10,10])
Objeto_FEM.perfilY(0.472,0.00098,0.00128,0,50)
Objeto_FEM.perfilX(0.125,0.0011,0.00114,0,1)
