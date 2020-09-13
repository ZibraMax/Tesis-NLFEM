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

class NoLocal(FEM):
    def __init__(this, geometria):
        super().__init__(2 * len(geometria.gdls))
        this.importarGeometria(geometria)

    def generarElementos(this, gauss=4):
        for d, t in zip(this.geometria.diccionarios, this.geometria.tipos):
            coords = np.array(this.geometria.gdls)[np.ix_(d)]
            if t == 'T1V' or t == 'T2V':
                if len(d) == 3:
                    k = np.array(d)
                    k1 = k * 2
                    k2 = k * 2 + 1
                    d1 = np.array(k1.tolist() + k2.tolist())
                    this.elementos.append(TriangularL(coords, d1, gauss=gauss))
                elif len(d) == 6:
                    k = np.array(d)
                    k1 = k * 2
                    k2 = k * 2 + 1
                    d1 = np.array(k1.tolist() + k2.tolist())
                    this.elementos.append(TriangularC(coords, d1, gauss=gauss))
                else:
                    raise 'algo anda mal'
            elif t == 'C1V':
                if len(d) == 8:
                    k = np.array(d)
                    k1 = k * 2
                    k2 = k * 2 + 1
                    d1 = np.array(k1.tolist() + k2.tolist())

                    this.elementos.append(SerendipityC(coords, np.array(d1), gauss=gauss))
                elif len(d) == 4:
                    k = np.array(d)
                    k1 = k * 2
                    k2 = k * 2 + 1
                    d1 = np.array(k1.tolist() + k2.tolist())
                    this.elementos.append(CuadrilateroL(coords, d1, gauss=gauss))
            else:
                raise Exception('No se pudo crear el elemento con coordenadas ' + format(coords))

    def calcularMatrices(this, E, v, Fx, Fy):
        distancia = lambda x0, y0, x1, y1: np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        count = 0
        for e in this.elementos:
            count += 1
            progressbar(count, len(this.elementos), prefix="Integrando elementos ", size=50)
            psis = e.psis
            dzpsis = e.dzpsis
            dnpsis = e.dnpsis
            gauss = e.gauss
            _J = e._J
            J = e.J
            x = e.Tx
            y = e.Ty

            n = len(psis(0, 0))

            Kuu = np.zeros([n, n])
            Kuv = np.zeros([n, n])
            Kvu = np.zeros([n, n])
            Kvv = np.zeros([n, n])

            Qu = np.zeros([n, 1])
            Qv = np.zeros([n, 1])

            Fu = np.zeros([n, 1])
            Fv = np.zeros([n, 1])

            K = np.zeros([2 * n, 2 * n])
            Q = np.zeros([2 * n, 1])
            F = np.zeros([2 * n, 1])
            C11 = E / (1 - v ** 2)
            C12 = v * E / (1 - v ** 2)
            C66 = E / 2 / (1 + v)
            for i in range(n):
                for j in range(n):
                    dfdx = lambda z, n, k: dzpsis(z, n)[k][0] * _J(z, n)[0][0] + dnpsis(z, n)[k][0] * _J(z, n)[0][1]
                    dfdy = lambda z, n, k: dzpsis(z, n)[k][0] * _J(z, n)[1][0] + dnpsis(z, n)[k][0] * _J(z, n)[1][1]
                    psi = lambda z, n, k: psis(z, n)[k][0]

                    EKuu = lambda z, n: (C11 * dfdx(z, n, i) * dfdx(z, n, j) + C66 * dfdy(z, n, i) * dfdy(z, n,
                                                                                                          j)) * np.linalg.det(
                        J(z, n))
                    EKuv = lambda z, n: (C12 * dfdx(z, n, i) * dfdy(z, n, j) + C66 * dfdy(z, n, i) * dfdx(z, n,
                                                                                                          j)) * np.linalg.det(
                        J(z, n))
                    EKvu = lambda z, n: (C12 * dfdy(z, n, i) * dfdx(z, n, j) + C66 * dfdx(z, n, i) * dfdy(z, n,
                                                                                                          j)) * np.linalg.det(
                        J(z, n))
                    EKvv = lambda z, n: (C11 * dfdy(z, n, i) * dfdy(z, n, j) + C66 * dfdx(z, n, i) * dfdx(z, n,
                                                                                                          j)) * np.linalg.det(
                        J(z, n))

                    Kuu[i, j] = e.intGauss2D(gauss, EKuu)
                    Kuv[i, j] = e.intGauss2D(gauss, EKuv)
                    Kvu[i, j] = e.intGauss2D(gauss, EKvu)
                    Kvv[i, j] = e.intGauss2D(gauss, EKvv)

                EFu = lambda z, n: (Fx(x(z, n), y(z, n))) * np.linalg.det(J(z, n))
                EFv = lambda z, n: (Fy(x(z, n), y(z, n))) * np.linalg.det(J(z, n))

                Fu[i] = e.intGauss2D(gauss, EFu)
                Fv[i] = e.intGauss2D(gauss, EFv)

            MD = np.linspace(0, 2 * n - 1, 2 * n).reshape([2, n]).astype(int)
            K[np.ix_(MD[0], MD[0])] = Kuu
            K[np.ix_(MD[0], MD[1])] = Kuv
            K[np.ix_(MD[1], MD[0])] = Kvu
            K[np.ix_(MD[1], MD[1])] = Kvv

            F[np.ix_(MD[0])] = Fu
            F[np.ix_(MD[1])] = Fv
            Knl = []
            countnl = 0
            for enl in this.elementos:
                countnl+=1
                psisnl = e.psis
                dzpsisnl = enl.dzpsis
                dnpsisnl = enl.dnpsis
                gaussnl = enl.gauss
                _Jnl = enl._J
                Jnl = enl.J
                xnl = enl.Tx
                ynl = enl.Ty

                nnl = len(psisnl(0, 0))

                Kuunl = np.zeros([nnl, nnl])
                Kuvnl = np.zeros([nnl, nnl])
                Kvunl = np.zeros([nnl, nnl])
                Kvvnl = np.zeros([nnl, nnl])

                Knlmn = np.zeros([2 * nnl, 2 * nnl])

                C11 = E / (1 - v ** 2)

                C12 = v * E / (1 - v ** 2)

                C66 = E / 2 / (1 + v)

                def A(X, XP, l=0.1,t=0.5):
                    l0 = (1)/(2*np.pi*l**2*t)
                    return l0 * np.exp(
                        -distancia(x(X[0], X[1]), y(X[0], X[1]), xnl(XP[0], XP[1]), ynl(XP[0], XP[1])) / l)

                for i in range(nnl):
                    for j in range(nnl):
                        dfdxnl = lambda z, n, k: dzpsisnl(z, n)[k][0] * _Jnl(z, n)[0][0] + dnpsisnl(z, n)[k][0] * \
                                                 _Jnl(z, n)[0][1]
                        dfdynl = lambda z, n, k: dzpsisnl(z, n)[k][0] * _Jnl(z, n)[1][0] + dnpsisnl(z, n)[k][0] * \
                                                 _Jnl(z, n)[1][1]
                        psinl = lambda z, n, k: psisnl(z, n)[k][0]

                        EKuunl = lambda z, n, znl, nnl: A([z, n], [znl, nnl]) * (
                                    C11 * dfdx(z, n, i) * dfdxnl(znl, nnl, j) + C66 * dfdy(z, n, i) * dfdynl(znl, nnl,
                                                                                                             j)) * np.linalg.det(
                            J(z, n)) * np.linalg.det(Jnl(znl, nnl))
                        EKuvnl = lambda z, n, znl, nnl: A([z, n], [znl, nnl]) * (
                                    C12 * dfdx(z, n, i) * dfdynl(znl, nnl, j) + C66 * dfdy(z, n, i) * dfdxnl(znl, nnl,
                                                                                                             j)) * np.linalg.det(
                            J(z, n)) * np.linalg.det(Jnl(znl, nnl))
                        EKvunl = lambda z, n, znl, nnl: A([z, n], [znl, nnl]) * (
                                    C12 * dfdy(z, n, i) * dfdxnl(znl, nnl, j) + C66 * dfdx(z, n, i) * dfdynl(znl, nnl,
                                                                                                             j)) * np.linalg.det(
                            J(z, n)) * np.linalg.det(Jnl(znl, nnl))
                        EKvvnl = lambda z, n, znl, nnl: A([z, n], [znl, nnl]) * (
                                    C11 * dfdy(z, n, i) * dfdynl(znl, nnl, j) + C66 * dfdx(z, n, i) * dfdxnl(znl, nnl,
                                                                                                             j)) * np.linalg.det(
                            J(z, n)) * np.linalg.det(Jnl(znl, nnl))

                        Kuunl[i, j] = e.intGauss4DNL(gaussnl, EKuunl)
                        Kuvnl[i, j] = e.intGauss4DNL(gaussnl, EKuvnl)
                        Kvunl[i, j] = e.intGauss4DNL(gaussnl, EKvunl)
                        Kvvnl[i, j] = e.intGauss4DNL(gaussnl, EKvvnl)

                Knlmn[np.ix_(MD[0], MD[0])] = Kuunl
                Knlmn[np.ix_(MD[0], MD[1])] = Kuvnl
                Knlmn[np.ix_(MD[1], MD[0])] = Kvunl
                Knlmn[np.ix_(MD[1], MD[1])] = Kvvnl
                Knl.append(Knlmn)
                progressbar(countnl, len(this.elementos), prefix="Integrando elementos No Locales   ", size=50)

            e.determinarMatrices(K, F, Q)
            e.KNLS = Knl


    def graficarSolucion(this, figsize=[12, 7], cmap='magma', linewidth=2, markersize=2, mask=None, mult=10):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()
        count = 0
        X_0 = []
        Y_0 = []
        X_F = []
        Y_F = []
        for e in this.elementos:
            count += 1

            Ue = this.U[np.ix_(e.gdl)]

            coords = e.coords

            e.coordsNuevas = np.array(coords) + Ue.T.reshape([2, len(e.psis(0, 0))]).T * mult
            coordsNuevas = e.coordsNuevas.tolist()
            coords = np.array(coords)

            X = np.array(e.coords)[:, 0].tolist()
            X.append(e.coords[0][0])
            Y = np.array(e.coords)[:, 1].tolist()
            Y.append(e.coords[0][1])
            e._coordenadas = np.array([X, Y]).T

            X = np.array(e.coordsNuevas)[:, 0].tolist()
            X.append(e.coordsNuevas[0][0])
            Y = np.array(e.coordsNuevas)[:, 1].tolist()
            Y.append(e.coordsNuevas[0][1])
            e._coordsNuevas = np.array([X, Y]).T

            X_0 = e._coordenadas[:, 0].tolist()
            Y_0 = e._coordenadas[:, 1].tolist()

            X_F = e._coordsNuevas[:, 0].tolist()
            Y_F = e._coordsNuevas[:, 1].tolist()

            progressbar(count, len(this.elementos), prefix="Graficando Solucion ", size=50)
            ax.plot(X_0, Y_0, 'o--', color='gray', zorder=0)
            ax.plot(X_F, Y_F, 'o-', color='black', zorder=3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Solución de elementos finitos')
        ax.legend(['Situacion Inicial', 'Deformada'])

    def defUnitariaX(this):
        count=0
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        xtotal = []
        ytotal = []
        ztotal = []
        markersize = 2
        cmap = 'magma'
        linewidth = 2
        for e in this.elementos:
            count += 1
            progressbar(count, len(this.elementos), prefix="Graficando derivada X ", size=50)
            x, y, u = postProcesoX(e,this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
            X = np.array(e._coordenadas)[:, 0]
            Y = np.array(e._coordenadas)[:, 1]
            ax.plot(X, Y, 'ko-', markersize=markersize, color='black', linewidth=linewidth)
        surf = ax.plot_trisurf(xtotal, ytotal, ztotal, cmap=cmap, zorder=1)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(r'$\frac{\partial U}{\partial X}$')

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
    this.U = lambda z, n: ex(this,z, n)
    for z, n in zip(Z, N):
        x.append(this.Tx(z, n)[0])
        y.append(this.Ty(z, n)[0])
        print(this.U(z, n))
        u.append(this.U(z, n))
    return x, y, u
def ex(e,z,n):
    dz = e.dzpsis(z, n)
    dn = e.dnpsis(z, n)
    result = []
    for i in range(len(dz)):
        result.append(e._J(z, n) @ np.array([[dz[i][0]], [dn[i][0]]]))
    result = np.array(result)
    dx = (e.Ue[np.ix_([0,2,4])].T @ result[:, 0])[0][0]
    return dx
