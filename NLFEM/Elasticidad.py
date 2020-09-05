import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import tri
import matplotlib as mpl
mpl.style.use('default')

from .FEM import *
from .FEM import progressbar as progressbar

from .SerendipityC import *
from .CuadrilateroL import *
from .TriangularC import *
from .TriangularL import *

class Elasticidad(FEM):
    def __init__(this,geometria):
        super().__init__(2*len(geometria.gdls))
        this.importarGeometria(geometria)
    
    def generarElementos(this,gauss=4):
        for d,t in zip(this.geometria.diccionarios,this.geometria.tipos):
            coords = np.array(this.geometria.gdls)[np.ix_(d)]
            if t == 'T1V':
                if len(d) == 3:
                    k = np.array(d)
                    k1 = k*2
                    k2 = k*2+1
                    d1 = np.array(k1.tolist()+k2.tolist())
                    this.elementos.append(TriangularL(coords,d1,gauss=gauss))
                elif len(d) == 6:
                    k = np.array(d)
                    k1 = k*2
                    k2 = k*2+1
                    d1 = np.array(k1.tolist()+k2.tolist())
                    this.elementos.append(TriangularC(coords,d1,gauss=gauss))
                else:
                    raise 'algo anda mal'
            elif t == 'C1V':
                if len(d) == 8:
                    k = np.array(d)
                    k1 = k*2
                    k2 = k*2+1
                    d1 = np.array(k1.tolist()+k2.tolist())
                    
                    this.elementos.append(SerendipityC(coords,np.array(d1),gauss=gauss))
                elif len(d) == 4:
                    k = np.array(d)
                    k1 = k*2
                    k2 = k*2+1
                    d1 = np.array(k1.tolist()+k2.tolist())
                    this.elementos.append(CuadrilateroL(coords,d1,gauss=gauss))
            else:
                raise Exception('No se pudo crear el elemento con coordenadas ' + format(coords))
    def calcularMatrices(this,E,v,Fx,Fy):
        count = 0
        for e in this.elementos:
            count += 1
            progressbar(count,len(this.elementos), prefix="Integrando elementos ", size=50)
            psis = e.psis
            dzpsis = e.dzpsis
            dnpsis = e.dnpsis
            gauss = e.gauss
            _J = e._J
            J = e.J
            x = e.Tx
            y = e.Ty
            
            n = len(psis(0,0))
            
            Kuu = np.zeros([n,n])
            Kuv = np.zeros([n,n])
            Kvu = np.zeros([n,n])
            Kvv = np.zeros([n,n])

            Qu = np.zeros([n,1])
            Qv = np.zeros([n,1])

            Fu = np.zeros([n,1])
            Fv = np.zeros([n,1])
            
            K = np.zeros([2*n,2*n])
            Q = np.zeros([2*n,1])
            F = np.zeros([2*n,1])
            C11 = E/(1-v**2)
            C12 = v*E/(1-v**2)
            C66 = E/2/(1+v)
            for i in range(n):
                for j in range(n):
                    dfdx = lambda z,n,k: dzpsis(z,n)[k][0]*_J(z,n)[0][0]+dnpsis(z,n)[k][0]*_J(z,n)[0][1]
                    dfdy = lambda z,n,k: dzpsis(z,n)[k][0]*_J(z,n)[1][0]+dnpsis(z,n)[k][0]*_J(z,n)[1][1]
                    psi = lambda z,n,k: psis(z,n)[k][0]
                    
                    EKuu = lambda z,n: (C11*dfdx(z,n,i)*dfdx(z,n,j)+C66*dfdy(z,n,i)*dfdy(z,n,j))*np.linalg.det(J(z,n))
                    EKuv = lambda z,n: (C12*dfdx(z,n,i)*dfdy(z,n,j)+C66*dfdy(z,n,i)*dfdx(z,n,j))*np.linalg.det(J(z,n))
                    EKvu = lambda z,n: (C12*dfdy(z,n,i)*dfdx(z,n,j)+C66*dfdx(z,n,i)*dfdy(z,n,j))*np.linalg.det(J(z,n))
                    EKvv = lambda z,n: (C11*dfdy(z,n,i)*dfdy(z,n,j)+C66*dfdx(z,n,i)*dfdx(z,n,j))*np.linalg.det(J(z,n))
                    
                    Kuu[i,j] = e.intGauss2D(gauss,EKuu)
                    Kuv[i,j] = e.intGauss2D(gauss,EKuv)
                    Kvu[i,j] = e.intGauss2D(gauss,EKvu)
                    Kvv[i,j] = e.intGauss2D(gauss,EKvv)
                    
                EFu = lambda z,n: (Fx(x(z,n),y(z,n)))*np.linalg.det(J(z,n))
                EFv = lambda z,n: (Fy(x(z,n),y(z,n)))*np.linalg.det(J(z,n))

                Fu[i] = e.intGauss2D(gauss,EFu)
                Fv[i] = e.intGauss2D(gauss,EFv)
            
            MD = np.linspace(0,2*n-1,2*n).reshape([2,n]).astype(int)
            K[np.ix_(MD[0],MD[0])] = Kuu
            K[np.ix_(MD[0],MD[1])] = Kuv
            K[np.ix_(MD[1],MD[0])] = Kvu
            K[np.ix_(MD[1],MD[1])] = Kvv
            
            F[np.ix_(MD[0])] = Fu
            F[np.ix_(MD[1])] = Fv
            e.determinarMatrices(K,F,Q)
    def graficarSolucion(this,figsize=[12,7],cmap='magma',linewidth=2,markersize=2,mask=None,mult=10):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        count = 0
        X_0 = []
        Y_0 = []
        X_F = []
        Y_F = []
        for e in this.elementos:
            count+=1
            
            Ue = this.U[np.ix_(e.gdl)]
            
            coords = e.coords
            
            e.coordsNuevas = np.array(coords) + Ue.T.reshape([2,len(e.psis(0,0))]).T*mult
            coordsNuevas = e.coordsNuevas.tolist()
            coords = np.array(coords)
            
            X = np.array(e.coords)[:,0].tolist()
            X.append(e.coords[0][0])
            Y = np.array(e.coords)[:,1].tolist()
            Y.append(e.coords[0][1])
            e._coordenadas = np.array([X,Y]).T
            
            X = np.array(e.coordsNuevas)[:,0].tolist()
            X.append(e.coordsNuevas[0][0])
            Y = np.array(e.coordsNuevas)[:,1].tolist()
            Y.append(e.coordsNuevas[0][1])
            e._coordsNuevas = np.array([X,Y]).T
            
            X_0 = e._coordenadas[:,0].tolist()
            Y_0 = e._coordenadas[:,1].tolist()
            
            X_F = e._coordsNuevas[:,0].tolist()
            Y_F = e._coordsNuevas[:,1].tolist()
            
            progressbar(count,len(this.elementos), prefix="Graficando Solucion ", size=50)
            ax.plot(X_0,Y_0,'o--',color='gray',zorder=0)
            ax.plot(X_F,Y_F,'o-',color='black',zorder=3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Solución de elementos finitos')
        ax.legend(['Situacion Inicial','Deformada'])