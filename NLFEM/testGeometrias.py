from Mesh import Geometria
from Mesh import delaunay
import matplotlib.pyplot as plt
vertices = [-1]
diccionarios = [[0,1,2],[0,2,4,3]]
gdls = [[1,3],[4,4],[3,2.8],[2,1],[4.5,1.5]]
tipos = ['T1V','T1V']
segmentos = [[0,3],[0,1]]
g = Geometria(vertices, diccionarios, gdls, tipos, segmentos=None)
g.dibujarse()


a=1
GEOMETRIA = delaunay.Delaunay1V([[0,0],[a,0],[a,a],[0,a]],delaunay._strdelaunay(a=0.1,o=2), plot=True,texto=10,bolita=10)
plt.show()