import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss as roots
from mpl_toolkits import mplot3d
from matplotlib import tri
import matplotlib as mpl
from matplotlib.lines import Line2D

mpl.style.use('default')


class Elemento:
    def __init__(this, coords, gdl=None, gauss=4):
        this.coords = coords
        this.__n = len(coords)
        this.gauss = gauss
        this.gdl = gdl
        this.transfCoords()
    def transfCoords(this):
        psis = this.psis
        dzpsis = this.dzpsis
        dnpsis = this.dnpsis
        coords = this.coords
        this.coordsNumpy = np.array(coords)
        dxdz = lambda z, n: this.coordsNumpy[:, 0] @ dzpsis(z, n)
        dydz = lambda z, n: this.coordsNumpy[:, 1] @ dzpsis(z, n)
        dxdn = lambda z, n: this.coordsNumpy[:, 0] @ dnpsis(z, n)
        dydn = lambda z, n: this.coordsNumpy[:, 1] @ dnpsis(z, n)
        this.Tx = lambda z, n: np.array(coords)[:, 0] @ psis(z, n)
        this.Ty = lambda z, n: np.array(coords)[:, 1] @ psis(z, n)
        this.J = lambda z, n: np.array([[dxdz(z, n)[0], dydz(z, n)[0]], [dxdn(z, n)[0], dydn(z, n)[0]]])
        this._J = lambda z, n: np.linalg.inv(this.J(z, n))

    def dibujarPsis(this):

        fig = plt.figure(figsize=[8.5, 5])

        ax = fig.add_subplot(projection='3d')

        Z = this._dominioNaturalZ
        N = this._dominioNaturalN
        psi = lambda z, n, k: this.psis(z, n)[k][0]
        count = 0
        pT = np.zeros([len(Z)])
        l = []
        for i in range(len(this.psis(0, 0))):
            x = []
            y = []
            p = []
            count += 1
            for z, n in zip(Z, N):
                x.append(this.Tx(z, n)[0])
                y.append(this.Ty(z, n)[0])
                p.append(psi(z, n, count - 1))
            pT += p
            surf = ax.plot_trisurf(x, y, p)
            l.append(r'$\Psi_{' + format(count) + '}$')
            surf._facecolors2d = surf._facecolors3d
            surf._edgecolors2d = surf._edgecolors3d
        l.append(r'$\sum_{i=1}^{' + format(count) + '}{\Psi_i}$')
        surf = ax.plot_trisurf(x, y, pT, alpha=0.4)
        surf._facecolors2d = surf._facecolors3d
        surf._edgecolors2d = surf._edgecolors3d
        ax.set_xlabel(r'$\zeta$')
        ax.set_ylabel(r'$\eta$')
        ax.set_title('Funciones de interpolación')
        ax.legend(l)
        return np.max(pT), np.min(pT), np.average(pT)

    def determinarMatrices(this, Ke, Fe, Qe):
        this.Ke = Ke
        this.Fe = Fe
        this.Qe = Qe

    def dibujarse(this):
        X = np.array(this.coords)[:, 0].tolist()
        X.append(this.coords[0][0])
        Y = np.array(this.coords)[:, 1].tolist()
        Y.append(this.coords[0][1])
        this._coordenadas = np.array([X, Y]).T
        fig = plt.figure()
        ax = fig.add_subplot()
        Z = np.zeros(len(X))
        surf = ax.plot(X, Y, '-o', color='black', markerSize=10)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')

    def graficaJacobiano(this, figsize=None):
        J = this.J
        x = this._dominioNaturalZ
        y = this._dominioNaturalN
        z = []
        for i, j in zip(x, y):
            z.append(np.linalg.det(J(i, j)))
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(projection='3d')
            surf = ax.plot_trisurf(x, y, np.array(z), cmap='magma')
            cbar = fig.colorbar(surf)
        ax.set_xlabel(r'$\zeta$')
        ax.set_ylabel(r'$\eta$')
        ax.set_zlabel(r'$|\mathcal{J}|$')
        ax.set_title('Determinante del Jacobiano de transformación')

    def darSolucion(this, U, graficar=False):
        X = np.array(this.coords)[:, 0].tolist()
        X.append(this.coords[0][0])
        Y = np.array(this.coords)[:, 1].tolist()
        Y.append(this.coords[0][1])
        this._coordenadas = np.array([X, Y]).T
        this.Ue = U[np.ix_(this.gdl)]
        this._Ue = this.Ue.T[0].tolist()
        this._Ue.append(this.Ue[0][0])
        Z = this._dominioNaturalZ
        N = this._dominioNaturalN
        x = []
        y = []
        u = []
        this.U = lambda z, n: this.Ue.T[0] @ this.psis(z, n)
        for z, n in zip(Z, N):
            x.append(this.Tx(z, n)[0])
            y.append(this.Ty(z, n)[0])
            u.append(this.U(z, n)[0])
        if graficar:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            surf = ax.plot_trisurf(x, y, u, cmap='magma')
            surf._facecolors2d = surf._facecolors3d
            surf._edgecolors2d = surf._edgecolors3d
            cbar = fig.colorbar(surf)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('u')
            ax.set_title('Solucion interpolada en el elemento')
        return x, y, u

    def _darSolucion(this, U, graficar=False):
        X = np.array(this.coords)[:, 0].tolist()
        X.append(this.coords[0][0])
        Y = np.array(this.coords)[:, 1].tolist()
        Y.append(this.coords[0][1])
        this._coordenadas = np.array([X, Y]).T
        Ue = U[np.ix_(this.gdl)]
        _Ue = Ue.T[0].tolist()
        _Ue.append(Ue[0][0])
        Z = this._dominioNaturalZ
        N = this._dominioNaturalN
        x = []
        y = []
        u = []
        U = lambda z, n: Ue.T[0] @ this.psis(z, n)
        for z, n in zip(Z, N):
            x.append(this.Tx(z, n)[0])
            y.append(this.Ty(z, n)[0])
            u.append(U(z, n)[0])
        if graficar:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            surf = ax.plot_trisurf(x, y, u, cmap='magma')
            surf._facecolors2d = surf._facecolors3d
            surf._edgecolors2d = surf._edgecolors3d
            cbar = fig.colorbar(surf)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('u')
            ax.set_title('Solucion interpolada en el elemento')
        return x, y, u

    def grad(this, z, n):
        dz = this.dzpsis(z, n)
        dn = this.dnpsis(z, n)
        result = []
        for i in range(len(dz)):
            result.append(this._J(z, n) @ np.array([[dz[i][0]], [dn[i][0]]]))
        result = np.array(result)
        dx = (this.Ue.T @ result[:, 0])[0][0]
        dy = (this.Ue.T @ result[:, 1])[0][0]
        return np.array([[dx], [dy]])

    def postProcesoX(this, U):

        X = np.array(this.coords)[:, 0].tolist()
        X.append(this.coords[0][0])
        Y = np.array(this.coords)[:, 1].tolist()
        Y.append(this.coords[0][1])
        this._coordenadas = np.array([X, Y]).T
        this.Ue = U[np.ix_(this.gdl)]
        this._Ue = this.Ue.T[0].tolist()
        this._Ue.append(this.Ue[0][0])
        Z = this._dominioNaturalZ
        N = this._dominioNaturalN
        x = []
        y = []
        u = []
        this.U = lambda z, n: this.grad(z, n)[0]
        for z, n in zip(Z, N):
            x.append(this.Tx(z, n)[0])
            y.append(this.Ty(z, n)[0])
            u.append(this.U(z, n)[0])
        return x, y, u

    def postProcesoY(this, U):

        X = np.array(this.coords)[:, 0].tolist()
        X.append(this.coords[0][0])
        Y = np.array(this.coords)[:, 1].tolist()
        Y.append(this.coords[0][1])
        this._coordenadas = np.array([X, Y]).T
        this.Ue = U[np.ix_(this.gdl)]
        this._Ue = this.Ue.T[0].tolist()
        this._Ue.append(this.Ue[0][0])
        Z = this._dominioNaturalZ
        N = this._dominioNaturalN
        x = []
        y = []
        u = []
        this.U = lambda z, n: this.grad(z, n)[1]
        for z, n in zip(Z, N):
            x.append(this.Tx(z, n)[0])
            y.append(this.Ty(z, n)[0])
            u.append(this.U(z, n)[0])
        return x, y, u

    def postProcesoGrad(this, U):
        X = np.array(this.coords)[:, 0].tolist()
        X.append(this.coords[0][0])
        Y = np.array(this.coords)[:, 1].tolist()
        Y.append(this.coords[0][1])
        this._coordenadas = np.array([X, Y]).T
        this.Ue = U[np.ix_(this.gdl)]
        this._Ue = this.Ue.T[0].tolist()
        this._Ue.append(this.Ue[0][0])
        Z = this._dominioNaturalZ
        N = this._dominioNaturalN
        x = []
        y = []
        u = []
        this.U = lambda z, n: np.sqrt((this.grad(z, n)[0][0]) ** 2 + (this.grad(z, n)[1][0]) ** 2)
        for z, n in zip(Z, N):
            x.append(this.Tx(z, n)[0])
            y.append(this.Ty(z, n)[0])
            u.append(this.U(z, n))
        return x, y, u

    def solucionEnPunto(this, x, y):
        """Interpola la solucion en un punto dado por parametro. El punto debe estar en el dominio del elemento"""
