from NLFEM.FEM1V import *
import NLFEM.Mesh as Mesh
from NLFEM.Mesh import delaunay

a = 0.3
b = 0.3
tw = 0.05
tf = 0.05

E = 200000
v = 0.27
G = E/(2*(1+v))

params = Mesh.delaunay._strdelaunay(constrained=True,delaunay=True,a='0.0005',o=2)
vertices = [[0,0],[a,0],[a,tf],[a/2+tw/2,tf],[a/2+tw/2,tf+b],[a,tf+b],[a,2*tf+b],[0,2*tf+b],[0,tf+b],[a/2-tw/2,tf+b],[a/2-tw/2,tf],[0,tf]]
geometria = Mesh.Delaunay1V(vertices, params, plot=True,texto=10,bolita=0)
plt.show()
plt.savefig('I_geometria.svg')
geometria.definirTodasCondiciones()
zanahorias = FEM1V(geometria)

zanahorias.generarElementos()

a11 = lambda x,y: 1
a12 = lambda x,y: 0
a21 = lambda x,y: 0
a22 = lambda x,y: 1
a00 = lambda x,y: 0

theta = 1
f = lambda x,y: 2*G*theta
zanahorias.definirCondicionesDeBorde(geometria.cbe)
zanahorias.solucionar(cmap='magma',markersize=1,linewidth=1,mask=vertices,a11=a11,a12=a12,a21=a21,a22=a22,a00=a00,f=f,plot=True)
plt.show()
plt.savefig('I.svg')
a=a