import NLFEM
from NLFEM.Mesh import Rect
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch,Ellipse
from matplotlib.lines import Line2D

legend_elements = [Patch(facecolor='orange',fill= False, hatch='\\\\' , edgecolor='r',
                         label='Centroide'),Line2D([0], [0], color='red', lw=4, label='Distancia Extra')
                   ]

a = 5
l = 0.1
GEOMETRIA = Rect.Rect("enmallado.txt")
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,0],[a,0])
GEOMETRIA.generarSegmentosDesdeCoordenadas([a,0],[a,a])
GEOMETRIA.generarSegmentosDesdeCoordenadas([a,a],[0,a])
GEOMETRIA.generarSegmentosDesdeCoordenadas([0,a],[0,0])
GEOMETRIA.dibujarse()

elemento = 63

base = GEOMETRIA.Centroides[elemento]

		


for i,ele in enumerate(GEOMETRIA.diccionariosnl[elemento]):
	e = GEOMETRIA._diccionarios[ele]
	coords = np.array(GEOMETRIA.gdls)[np.ix_(e)]
	coords = np.array(coords.tolist() + [coords[0].tolist()])
	X = coords[:, 0]
	Y = coords[:, 1]
	plt.fill(X, Y, color='red',zorder=-100,alpha=0.6)
for i,e in enumerate(GEOMETRIA._diccionarios):
	r = ((GEOMETRIA.Centroides[i][0]-base[0])**2+(GEOMETRIA.Centroides[i][1]-base[1])**2)**0.5
	if r <= 6*l:
		coords = np.array(GEOMETRIA.gdls)[np.ix_(e)]
		coords = np.array(coords.tolist() + [coords[0].tolist()])
		X = coords[:, 0]
		Y = coords[:, 1]
		plt.fill(X, Y, color='orange',fill=False, hatch='\\\\',zorder=-100)
plt.gca().add_patch(Ellipse(base, 12*l, 12*l,fill= False, edgecolor='b'))
plt.legend(handles=legend_elements)
plt.show()