import NLFEM
from NLFEM.Mesh import Rect
import subprocess
import matplotlib.pyplot as plt
import numpy as np

RUTA_M = 'NLFEM_C++/MATRICESPolimorfismo'

E = 2.1*10**6
V = 0.2
u = 0.001
a = 5
t = 0.5
l = 0.1

GEOMETRIA = Rect.Rect("NLFEM_C++/input.txt")
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,0],[a,0])
GEOMETRIA.generarSegmentosDesdeCoordenadas([a,0],[a,a])
GEOMETRIA.generarSegmentosDesdeCoordenadas([a,a],[0,a])
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,a],[0,0])
condiciones_borde_escenciales = GEOMETRIA.generarCBdesdeBorde(3,[0,0])+GEOMETRIA.generarCBdesdeBordeX(1,u)
GEOMETRIA.dibujarse()
plt.show()

Objeto_FEM = NLFEM.NoLocal(GEOMETRIA)
Objeto_FEM.z1 = 0.5
Objeto_FEM.generarElementos()
Objeto_FEM.importarMatricesCarpeta(RUTA_M)
Objeto_FEM.definirCondicionesDeBorde(condiciones_borde_escenciales)
Objeto_FEM.ensamblar()
Objeto_FEM.condicionesFrontera()
Objeto_FEM.solucionarSistemaEcuaciones()

Objeto_FEM_Local = NLFEM.NoLocal(GEOMETRIA)
Objeto_FEM_Local.z1 = 1
Objeto_FEM_Local.generarElementos()
Objeto_FEM_Local.importarMatricesCarpeta(RUTA_M)
Objeto_FEM_Local.definirCondicionesDeBorde(condiciones_borde_escenciales)
Objeto_FEM_Local.ensamblar()
Objeto_FEM_Local.condicionesFrontera()
Objeto_FEM_Local.solucionarSistemaEcuaciones()

Objeto_FEM.defUnitariaX([10,10])
plt.show()

Objeto_FEM.perfilY(0.019,0.00016,0.0004,0,a,label='No Local')
Objeto_FEM_Local.perfilY(0.019,0.00016,0.0004,0,a,acum=True,label='Local')
plt.show()

Objeto_FEM.perfilY(2.519,0.00018,0.00028,0,a,label='No Local')
Objeto_FEM_Local.perfilY(2.519,0.00018,0.00028,0,a,acum=True,label='Local')
plt.show()

Objeto_FEM.perfilX(0.019,0.00016,0.0004,0,a,label='No Local')
Objeto_FEM_Local.perfilX(0.019,0.00016,0.0004,0,a,acum=True,label='Local')
plt.show()

Objeto_FEM.perfilX(2.519,0.00016,0.0004,0,a,label='No Local')
Objeto_FEM_Local.perfilX(2.519,0.00016,0.0004,0,a,acum=True,label='Local')
plt.show()
