import NLFEM
from NLFEM.Mesh import Rect
import matplotlib.pyplot as plt

RUTA_M = 'MATRICES'

E = 2.1*10**6
V = 0.2
u = 0.001
a = 5 #Base del rectángulo
t = 0.5 #Grosor
l = 0.1 

GEOMETRIA = Rect.Rect("enmallado.txt")
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,0],[a,0])
GEOMETRIA.generarSegmentosDesdeCoordenadas([a,0],[a,a])
GEOMETRIA.generarSegmentosDesdeCoordenadas([a,a],[0,a])
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,a],[0,0])
GEOMETRIA.dibujarse()

condiciones_borde_escenciales = GEOMETRIA.generarCBdesdeBorde(3,[0,0])+GEOMETRIA.generarCBdesdeBordeX(1,u)

Objeto_FEM = NLFEM.NoLocal(GEOMETRIA)
Objeto_FEM.z1 = 0.5
# Objeto_FEM.moduloIntegrador("enmallado.txt",3,E,V,t,l,RUTA_M,tfa=1)
Objeto_FEM.generarElementos()
Objeto_FEM.importarMatricesCarpeta(RUTA_M)
Objeto_FEM.definirCondicionesDeBorde(condiciones_borde_escenciales)
Objeto_FEM.ensamblar()
Objeto_FEM.condicionesFrontera()
Objeto_FEM.solucionarSistemaEcuaciones()


import matplotlib.pyplot as plt
Objeto_FEM.defUnitariaX([10,10])
plt.show()
Objeto_FEM.perfilY(0.019,0.00016,0.0004,0,a,acum=True,label='No Local')
plt.show()
