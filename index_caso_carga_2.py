import NLFEM
import numpy as np
import matplotlib.pyplot as plt
from NLFEM.Mesh import Rect
import os

import copy

MATRICES_DESDE = 'C++'
RUTA_M = ''
if MATRICES_DESDE == 'FORTRAN':
	RUTA_M = 'NLFEM_FORTRAN'
elif MATRICES_DESDE == 'C++':
	RUTA_M = 'NLFEM_C++/NLFEM l 0.1/MATRICES'
else:
	print('Arturo por favor deja de jugar')

def condicionCarga(x,U0,a):
    m = U0/(a/2)
    if x<=a/2:
        return x*m
    else:
        return -x*m+U0*2

def postProcesoX(this, U):
    this.Ue = U[np.ix_(this.gdl)]
    this._Ue = this.Ue.T[0].tolist()
    this._Ue.append(this.Ue[0][0])
    Z = this._dominioNaturalZ
    N = this._dominioNaturalN
    x = []
    y = []
    u = []
    this.U = lambda z, n: grad(this, z, n)[0]
    for z, n in zip(Z, N):
        x.append(this.Tx(z, n)[0])
        y.append(this.Ty(z, n)[0])
        u.append(this.U(z, n)[0])
    return x, y, u

def postProcesoXNodos(this, U):
    this.Ue = U[np.ix_(this.gdl)]
    this._Ue = this.Ue.T[0].tolist()
    this._Ue.append(this.Ue[0][0])
    Z = this.ZNatural
    N = this.NNatural
    u = []
    this.U = lambda z, n: grad(this, z, n)[0]
    for z, n in zip(Z, N):
        u.append(this.U(z, n)[0])
    return np.array(u).reshape([len(Z),1])

def postProcesoY(this, U):
    this.Ue = U[np.ix_(this.gdl)]
    this._Ue = this.Ue.T[0].tolist()
    this._Ue.append(this.Ue[0][0])
    Z = this._dominioNaturalZ
    N = this._dominioNaturalN
    x = []
    y = []
    u = []
    this.U = lambda z, n: grad(this, z, n)[1]
    for z, n in zip(Z, N):
        x.append(this.Tx(z, n)[0])
        y.append(this.Ty(z, n)[0])
        u.append(this.U(z, n)[0])
    return x, y, u

def grad(this, z, n):
    dz = this.dzpsis(z, n)
    dn = this.dnpsis(z, n)
    result = []
    for i in range(len(dz)):
        result.append(this._J(z, n) @ np.array([[dz[i][0]], [dn[i][0]]]))
    result = np.array(result)
    n = len(result)
    U = np.linspace(0,n-1,n).astype(int)
    V = U*2+1
    dx = (this.Ue[np.ix_(U)].T @ result[:, 0])[0]
    dy = (this.Ue[np.ix_(V)].T @ result[:, 1])[0]
    return np.array([dx, dy])
def defUnitariaX(this,figsize):
    count=0
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')
    markersize = 2
    cmap = 'magma'
    linewidth = 2
    xtotal = []
    ytotal = []
    ztotal = []
    count = 0
    for e in this.elementos:
        count+=1
        x,y,u = postProcesoX(e,this.U)
        xtotal.extend(x)
        ytotal.extend(y)
        ztotal.extend(u)
    surf = ax.plot_trisurf(xtotal, ytotal, ztotal,cmap=cmap,zorder=1)
    cbar = fig.colorbar(surf)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel(r'$\varepsilon x$')
    ax.set_title(r'$\frac{\partial U}{\partial X}$')
    return xtotal,ytotal,ztotal
def defUnitariaY(this,figsize):
    count=0
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection='3d')
    markersize = 2
    cmap = 'magma'
    linewidth = 2
    xtotal = []
    ytotal = []
    ztotal = []
    count = 0
    for e in this.elementos:
        count+=1
        x,y,u = postProcesoY(e,this.U)
        xtotal.extend(x)
        ytotal.extend(y)
        ztotal.extend(u)
    surf = ax.plot_trisurf(xtotal, ytotal, ztotal,cmap=cmap,zorder=1)
    cbar = fig.colorbar(surf)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel(r'$\varepsilon x$')
    ax.set_title(r'$\frac{\partial U}{\partial Y}$')

def darDeformacionXY(this, U, x, y, n=100):
    this.Ue = U[np.ix_(this.gdl)]
    this.U = lambda z, n: this.Ue.T[0] @ this.psis(z, n)
    zeta,eta = this.mappingInverso(x, y, n).reshape(1,2)[0].tolist()
    return this.U(zeta,eta)

PATH = os.getcwd()
nombre = "\\NLFEM\\Mesh\\input.txt"
GEOMETRIA = Rect.Rect(PATH+nombre,10,10)
GEOMETRIA.segmentos = [[0,60],[60,2820],[2820,2760],[2760,0]]
GEOMETRIA.cbe = []

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
def generarCBdesdeBordeXCaso2(this,borde,valor,tborde):
    cb = []
    nodos = this.darNodosCB(borde)
    cbe = np.zeros([len(nodos), 2])
    cbe[:, 0] = nodos*2
    cbe[:, 1] = [condicionCarga(x,valor,tborde) for x in np.array(GEOMETRIA.gdls)[np.ix_(nodos)][:,1]]
    cb += cbe.tolist()
    return cb
def generarCBdesdeBorde(this, borde, valor=[0,0]):
    return generarCBdesdeBordeX(this, borde, valor[0])+generarCBdesdeBordeY(this, borde, valor[1])

a = 5
E = 2.1*10**6
V = 0.2
u = 0.005
Objeto_FEM = NLFEM.NoLocal(GEOMETRIA)
GEOMETRIA.cbe = generarCBdesdeBorde(GEOMETRIA,3,[0,0])+generarCBdesdeBordeXCaso2(GEOMETRIA,1,u,a)
Objeto_FEM.definirCondicionesDeBorde(GEOMETRIA.cbe)
Objeto_FEM.generarElementos()
for i,e in enumerate(Objeto_FEM.elementos):
    n = len(e.coords)
    K = np.loadtxt(RUTA_M+'/Elemento'+format(i+1)+'/KL_'+format(i+1)+'.csv',delimiter=',').reshape([16,16])
    Q = np.zeros([2 * n, 1])
    F = np.zeros([2 * n, 1])
    e.determinarMatrices(K,F,Q)
    knls = []
    j=0
    KNLS = np.loadtxt(RUTA_M+'/Elemento'+format(i+1)+'/KNLS.csv',delimiter=',')
    for _ in e.elementosnl:
        KNL = KNLS[j].reshape([16,16]).T
        j+=1
        knls.append(KNL)
    e.KNLS = knls
    print('Cargando matrices del elemento ' + format(i))
PUNTOS = [-np.sqrt(3.0/5.0),0,np.sqrt(3.0/5.0)]
for this in Objeto_FEM.elementos:
    this._dominioNaturalZ = PUNTOS
    this._dominioNaturalN = PUNTOS
z1=0.5
Objeto_FEM.z1 = z1
Objeto_FEM.K = np.zeros([Objeto_FEM.n,Objeto_FEM.n])
Objeto_FEM.F = np.zeros([Objeto_FEM.n,1])
Objeto_FEM.Q = np.zeros([Objeto_FEM.n,1])
Objeto_FEM.U = np.zeros([Objeto_FEM.n,1])
Objeto_FEM.S = np.zeros([Objeto_FEM.n,1])
Objeto_FEM.definirCondicionesDeBorde(Objeto_FEM.geometria.cbe)
Objeto_FEM.ensamblar()
Objeto_FEM.condicionesFrontera(Objeto_FEM.cbe,Objeto_FEM.cbn)
Objeto_FEM.solucionarSistemaEcuaciones()

Objeto_FEMLOCAL = copy.deepcopy(Objeto_FEM)
Objeto_FEMLOCAL.z1 = 1
Objeto_FEMLOCAL.K = np.zeros([Objeto_FEMLOCAL.n,Objeto_FEMLOCAL.n])
Objeto_FEMLOCAL.F = np.zeros([Objeto_FEMLOCAL.n,1])
Objeto_FEMLOCAL.Q = np.zeros([Objeto_FEMLOCAL.n,1])
Objeto_FEMLOCAL.U = np.zeros([Objeto_FEMLOCAL.n,1])
Objeto_FEMLOCAL.S = np.zeros([Objeto_FEMLOCAL.n,1])
Objeto_FEMLOCAL.definirCondicionesDeBorde(Objeto_FEMLOCAL.geometria.cbe)
Objeto_FEMLOCAL.ensamblar()
Objeto_FEMLOCAL.condicionesFrontera(Objeto_FEMLOCAL.cbe,Objeto_FEMLOCAL.cbn)
Objeto_FEMLOCAL.solucionarSistemaEcuaciones()

# [XNoLocal, YNoLocal, ZNoLocal] = defUnitariaX(Objeto_FEM,[10,10])
# plt.show()

def perfilX(x,yi=0.00016,yf=0.0004):
    _Y = np.linspace(0,5,1000).tolist()
    Z = []
    Y = []
    Zlocal = []
    for y in _Y:
        for i in range(len(Objeto_FEM.elementos)):
            e = Objeto_FEM.elementos[i]
            elocal = Objeto_FEMLOCAL.elementos[i]
            if e.estaDentro(x,y):
                U = postProcesoXNodos(e,Objeto_FEM.U)
                Y.append(y)
                Z.append(e.darSolucionXY(U, x, y, n=100))
                U = postProcesoXNodos(elocal,Objeto_FEMLOCAL.U)
                Zlocal.append(elocal.darSolucionXY(U, x, y, n=100))
                break
    plt.plot(Y,Z)
    plt.plot(Y,Zlocal)

    plt.xlim(0, 5)
    plt.ylim(yi, yf)

    plt.grid()
    plt.title('Perfil a x='+format(x))
    plt.xlabel('y')
    plt.ylabel(r'$\frac{\partial U}{\partial x}=\varepsilon_x$')
    plt.legend(['No Local','Local'])
    plt.savefig('Caso2deformacionesXEnX'+format(x)+'.png',transparent=True)
    plt.show()
def perfilY(y,yi=0.00016,yf=0.0004):
    _X = np.linspace(0,5,1000).tolist()
    Z = []
    X = []
    Zlocal = []
    for x in _X:
        for i in range(len(Objeto_FEM.elementos)):
            e = Objeto_FEM.elementos[i]
            elocal = Objeto_FEMLOCAL.elementos[i]
            if e.estaDentro(x,y):
                U = postProcesoXNodos(e,Objeto_FEM.U)
                X.append(x)
                Z.append(e.darSolucionXY(U, x, y, n=100))
                U = postProcesoXNodos(elocal,Objeto_FEMLOCAL.U)
                Zlocal.append(elocal.darSolucionXY(U, x, y, n=100))
                break
    plt.plot(X,Z)
    plt.plot(X,Zlocal)

    plt.xlim(0, 5)
    plt.ylim(yi, yf)

    plt.grid()
    plt.title('Perfil a y='+format(y))
    plt.xlabel('x')
    plt.ylabel(r'$\frac{\partial U}{\partial x}=\varepsilon_x$')
    plt.legend(['No Local','Local'])
    plt.savefig('Caso2deformacionesXEnY'+format(y)+'.png',transparent=True)
    plt.show()
#Vamos a ver, perfiles:

#Figura 9
perfilY(0.250,-0.0008,0.0008)
perfilY(1.519,0.0004,0.0008)
perfilY(2.519,0,0.004)

#Figura 10
perfilX(0.019,0.0004,0.0012)
perfilX(2.519,0,0.001)
perfilX(4.981,-0.001,0.004)
