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

class FEM1V(FEM):
    def __init__(this,geometria):
        super().__init__(len(geometria.gdls))
        this.importarGeometria(geometria)
        this.cbe = geometria.cbe
        this.cbn = geometria.cbn
        
    def generarElementos(this,gauss=4):
        for d,t in zip(this.geometria.diccionarios,this.geometria.tipos):
            coords = np.array(this.geometria.gdls)[np.ix_(d)]
            if t == 'T1V' or t == 'T2V':
                if len(d) == 3:
                    this.elementos.append(TriangularL(coords,np.array(d),gauss=gauss))
                elif len(d) == 6:
                    this.elementos.append(TriangularC(coords,np.array(d),gauss=gauss))
                else:
                    raise 'algo anda mal'
            elif t == 'C1V':
                if len(d) == 8:
                    this.elementos.append(SerendipityC(coords,np.array(d),gauss=gauss))
                elif len(d) == 4:
                    this.elementos.append(CuadrilateroL(coords,np.array(d),gauss=gauss))
    def calcularMatrices(this,a11,a12,a21,a22,a00,f):
        count = 0
        for e in this.elementos:
            count += 1
            progressbar(count,len(this.elementos), prefix="Integrando elementos ", size=50)
            psis = e.psis
            dzpsi = e.dzpsi
            dnpsi = e.dnpsi
            gauss = e.gauss
            _J = e._J
            J = e.J
            x = e.Tx
            y = e.Ty
            n = len(psis(0,0))
            K = np.zeros([n,n])
            Q = np.zeros([n,1])
            F = np.zeros([n,1])
            psi = lambda z,n,k: e.psi[k](z,n)
            for i in range(n):
                for j in range(n):
                    def EK(z,n):
                        X = x(z,n)
                        Y = y(z,n)
                        jacobiano = _J(z,n)
                        detjac = np.linalg.det(J(z,n))
                        dfdx_i = dzpsi[i](n, z) * jacobiano[0][0] + dnpsi[i](n, z) * jacobiano[0][1]
                        dfdy_i = dzpsi[i](n, z) * jacobiano[1][0] + dnpsi[i](n, z) * jacobiano[1][1]

                        dfdx_j = dzpsi[j](n, z) * jacobiano[0][0] + dnpsi[j](n, z) * jacobiano[0][1]
                        dfdy_j = dzpsi[j](n, z) * jacobiano[1][0] + dnpsi[j](n, z) * jacobiano[1][1]

                        psi_i = e.psi[i](z,n)
                        psi_j = e.psi[j](z,n)
                        return (a11(X,Y)*dfdx_i*dfdx_j + a12(X,Y)*dfdx_i*dfdy_j + a21(X,Y)*dfdy_i*dfdx_j + a22(X,Y)*dfdy_i*dfdy_j +a00(X,Y)*psi_i*psi_j )*detjac
                    K[i,j] = e.intGauss2D(gauss,EK)
                EF = lambda z,n: f(x(z,n),y(z,n))*psi(z,n,i)*np.linalg.det(J(z,n))
                F[i] = e.intGauss2D(gauss,EF)
            e.determinarMatrices(K,F,Q)
    def graficarSolucion(this,figsize=[17,10],cmap='magma',linewidth=2,markersize=2,mask=None):
        xtotal = []
        ytotal = []
        ztotal = []
        fig = plt.figure(figsize=figsize,constrained_layout=True)
        gs = fig.add_gridspec(2, 2)
        ax = fig.add_subplot(gs[0, 0],projection='3d')
        count = 0
        for e in this.elementos:
            count+=1
            x,y,u = e.darSolucion(this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
            X = np.array(e._coordenadas)[:,0]
            Y = np.array(e._coordenadas)[:,1]
            U = e._Ue
            ax.plot(X,Y,U, 'ko-',markersize=markersize,color='black',linewidth=linewidth)
            progressbar(count,len(this.elementos), prefix="Graficando Solucion ", size=50)
        surf = ax.plot_trisurf(xtotal, ytotal, ztotal,cmap=cmap)
        surf._facecolors2d=surf._facecolors3d
        surf._edgecolors2d=surf._edgecolors3d
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u')
        ax.set_title('Solucion de elementos finitos')
        
        ax = fig.add_subplot(gs[0, 1])

        count = 0
        for e in this.elementos:
            progressbar(count,len(this.elementos), prefix="Graficando contornos ", size=50)
            count+=1
            X = np.array(e._coordenadas)[:,0]
            Y = np.array(e._coordenadas)[:,1]
            ax.plot(X,Y,'ko-',markersize=markersize,color='black',linewidth=linewidth)
        surf = ax.tricontourf(xtotal, ytotal, ztotal,cmap=cmap,zorder=1)
        if not mask == None:
            cornersnt = np.array(mask[::-1])

            xmin = np.min(cornersnt[:,0])
            xmax = np.max(cornersnt[:,0])

            ymin = np.min(cornersnt[:,1])
            ymax = np.max(cornersnt[:,1])

            Xs = [xmin,xmax,xmax,xmin]+cornersnt[:,0].tolist()
            Ys = [ymin,ymin,ymax,ymax]+cornersnt[:,1].tolist()
            ax.fill(Xs,Ys,color='white',zorder=30)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axes.set_aspect('equal')
        ax.set_title('Solucion de elementos finitos')
        
        ax = fig.add_subplot(gs[1, 0])
        xtotal = []
        ytotal = []
        ztotal = []
        count = 0
        for e in this.elementos:
            count+=1
            progressbar(count,len(this.elementos), prefix="Graficando derivada x ", size=50)
            x,y,u = e.postProcesoX(this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
            X = np.array(e._coordenadas)[:,0]
            Y = np.array(e._coordenadas)[:,1]
            ax.plot(X,Y,'ko-',markersize=markersize,color='black',linewidth=linewidth)
        surf = ax.tricontourf(xtotal, ytotal, ztotal,cmap=cmap,zorder=1)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axes.set_aspect('equal')
        ax.set_title(r'$\frac{\partial U}{\partial X}$')
        if not mask== None:
            cornersnt = np.array(mask[::-1])

            xmin = np.min(cornersnt[:,0])
            xmax = np.max(cornersnt[:,0])

            ymin = np.min(cornersnt[:,1])
            ymax = np.max(cornersnt[:,1])

            Xs = [xmin,xmax,xmax,xmin]+cornersnt[:,0].tolist()
            Ys = [ymin,ymin,ymax,ymax]+cornersnt[:,1].tolist()
            ax.fill(Xs,Ys,color='white',zorder=30)
        
        ax = fig.add_subplot(gs[1, 1])
        xtotal = []
        ytotal = []
        ztotal = []
        count=0
        for e in this.elementos:
            count+=1
            progressbar(count,len(this.elementos), prefix="Graficando derivada y ", size=50)
            x,y,u = e.postProcesoY(this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
            X = np.array(e._coordenadas)[:,0]
            Y = np.array(e._coordenadas)[:,1]
            ax.plot(X,Y,'ko-',markersize=markersize,color='black',linewidth=linewidth)
        surf = ax.tricontourf(xtotal, ytotal, ztotal,cmap=cmap,zorder=1)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axes.set_aspect('equal')
        ax.set_title(r'$\frac{\partial U}{\partial Y}$')
        if not mask== None:
            cornersnt = np.array(mask[::-1])

            xmin = np.min(cornersnt[:,0])
            xmax = np.max(cornersnt[:,0])

            ymin = np.min(cornersnt[:,1])
            ymax = np.max(cornersnt[:,1])

            Xs = [xmin,xmax,xmax,xmin]+cornersnt[:,0].tolist()
            Ys = [ymin,ymin,ymax,ymax]+cornersnt[:,1].tolist()
            ax.fill(Xs,Ys,color='white',zorder=30)
    def graficarMagnitudGradiente(this,cmap='magma',figsize=None,mask=None):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        xtotal = []
        ytotal = []
        ztotal = []
        for e in this.elementos:
            x,y,u = e.postProcesoGrad(this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
        surf = ax.tricontourf(xtotal, ytotal, ztotal,cmap=cmap,zorder=2.5)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(r'$\sqrt{\frac{\partial U}{\partial X}^2+\frac{\partial U}{\partial X}^2}$')
        
        if not mask== None:
            cornersnt = np.array(mask[::-1])

            xmin = np.min(cornersnt[:,0])
            xmax = np.max(cornersnt[:,0])

            ymin = np.min(cornersnt[:,1])
            ymax = np.max(cornersnt[:,1])

            Xs = [xmin,xmax,xmax,xmin]+cornersnt[:,0].tolist()
            Ys = [ymin,ymin,ymax,ymax]+cornersnt[:,1].tolist()
            ax.fill(Xs,Ys,color='white',zorder=30)
    def graficarSolucionFast(this,U,figsize=[12,7],cmap='magma',linewidth=2,markersize=2,mask=None,name='Solucion de elementos finitos'):
        xtotal = []
        ytotal = []
        ztotal = []
        fig = plt.figure(figsize=figsize,constrained_layout=True)
        ax = fig.add_subplot()
        for e in this.elementos:
            x,y,u = e._darSolucion(U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
            X = np.array(e._coordenadas)[:,0]
            Y = np.array(e._coordenadas)[:,1]
            ax.plot(X,Y,'ko-',markersize=markersize,color='black',linewidth=linewidth)
        surf = ax.tricontourf(xtotal, ytotal, ztotal,cmap=cmap,zorder=1)
        if not mask == None:
            cornersnt = np.array(mask[::-1])

            xmin = np.min(cornersnt[:,0])
            xmax = np.max(cornersnt[:,0])

            ymin = np.min(cornersnt[:,1])
            ymax = np.max(cornersnt[:,1])

            Xs = [xmin,xmax,xmax,xmin]+cornersnt[:,0].tolist()
            Ys = [ymin,ymin,ymax,ymax]+cornersnt[:,1].tolist()
            ax.fill(Xs,Ys,color='white',zorder=30)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(name)