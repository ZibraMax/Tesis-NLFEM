from .TriangularL import *
from .TriangularC import *
from .CuadrilateroL import *
from .SerendipityC import *
from .FEM import progressbar as progressbar
from .FEM import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import tri
import matplotlib as mpl
import subprocess
mpl.style.use('default')


class NoLocal(FEM):
    def __init__(this, geometria):
        super().__init__(2 * len(geometria.gdls))
        this.importarGeometria(geometria)

    def generarElementos(this, gauss=3):
        i = -1
        for d, t in zip(this.geometria.diccionarios, this.geometria.tipos):
            i += 1
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
            elif t == 'C1V' or t == 'C2V':
                if len(d) == 8:
                    k = np.array(d)
                    k1 = k * 2
                    k2 = k * 2 + 1
                    d1 = np.array(k1.tolist() + k2.tolist())

                    this.elementos.append(SerendipityC(
                        coords, np.array(d1), gauss=gauss))
                elif len(d) == 4:
                    k = np.array(d)
                    k1 = k * 2
                    k2 = k * 2 + 1
                    d1 = np.array(k1.tolist() + k2.tolist())
                    this.elementos.append(
                        CuadrilateroL(coords, d1, gauss=gauss))
            else:
                raise Exception(
                    'No se pudo crear el elemento con coordenadas ' + format(coords))
            try:
                this.elementos[-1].elementosnl = this.geometria.diccionariosnl[i]
            except:
                this.elementos[-1].elementosnl = None
        if this.elementos[-1].elementosnl == None:
            for e in this.elementos:
                e.elementosnl = np.linspace(
                    0, len(this.elementos)-1, len(this.elementos)).astype(int).tolist()
        PUNTOS = [-np.sqrt(3.0/5.0), 0, np.sqrt(3.0/5.0)]
        for e in this.elementos:
            e._dominioNaturalZ = PUNTOS
            e._dominioNaturalN = PUNTOS
        this.K = np.zeros([this.n, this.n])
        this.F = np.zeros([this.n, 1])
        this.Q = np.zeros([this.n, 1])
        this.U = np.zeros([this.n, 1])
        this.S = np.zeros([this.n, 1])

    def importarMatricesCarpeta(this, RUTA_M):
        print('Importando matrices...')
        for i, e in tqdm(enumerate(this.elementos), total=len(this.elementos), unit='Elemento'):
            n = len(e.coords)
            K = np.loadtxt(RUTA_M+'/Elemento'+format(i+1)+'/KL_'+format(i+1) +
                           '.csv', delimiter=',').reshape([int(2*n), int(2*n)])
            Q = np.zeros([2 * n, 1])
            F = np.zeros([2 * n, 1])
            e.determinarMatrices(K, F, Q)
            knls = []
            j = 0
            KNLS = np.loadtxt(RUTA_M+'/Elemento'+format(i+1) +
                              '/KNLS.csv', delimiter=',')
            for enl in e.elementosnl:
                KNL = KNLS[j].reshape([int(2*n), int(2*n)]).T
                j += 1
                knls.append(KNL)
            e.KNLS = knls
        # PUNTOS = [-np.sqrt(3.0/5.0),0,np.sqrt(3.0/5.0)]
        # for e in this.elementos:
            # e._dominioNaturalZ = PUNTOS
            # e._dominioNaturalN = PUNTOS
        this.K = np.zeros([this.n, this.n])
        this.F = np.zeros([this.n, 1])
        this.Q = np.zeros([this.n, 1])
        this.U = np.zeros([this.n, 1])
        this.S = np.zeros([this.n, 1])
        print('Las matrices se han cargado con éxito')

    def calcularMatrices(this, E, v, Fx, Fy, t):
        def distancia(x0, y0, x1, y1): return np.sqrt(
            (x1 - x0) ** 2 + (y1 - y0) ** 2)
        count = 0
        for e in this.elementos:
            count += 1
            progressbar(count, len(this.elementos),
                        prefix="Integrando elementos ", size=50)

            def psi(z, n, k): return e.psi[k](z, n)
            dzpsi = e.dzpsi
            dnpsi = e.dnpsi
            gauss = e.gauss
            _J = e._J
            J = e.J
            x = e.Tx
            y = e.Ty

            n = len(e.coords)

            Kuu = np.zeros([n, n])
            Kuv = np.zeros([n, n])
            Kvu = np.zeros([n, n])
            Kvv = np.zeros([n, n])

            Fu = np.zeros([n, 1])
            Fv = np.zeros([n, 1])

            K = np.zeros([2 * n, 2 * n])
            Q = np.zeros([2 * n, 1])
            F = np.zeros([2 * n, 1])
            C11 = E / (1 - v ** 2)
            C12 = v * E / (1 - v ** 2)
            C66 = E / 2 / (1 + v)
            print(J(0, 0))
            for i in range(n):
                for j in range(n):
                    def dfdx(n, z, k):
                        jacobiano = _J(n, z)
                        return dzpsi[k](n, z) * jacobiano[0][0] + dnpsi[k](n, z) * jacobiano[0][1]

                    def dfdy(n, z, k):
                        jacobiano = _J(n, z)
                        return dzpsi[k](n, z) * jacobiano[1][0] + dnpsi[k](n, z) * jacobiano[1][1]

                    def EKuu(z, n): return (C11 * dfdx(z, n, i) * dfdx(z, n, j) + C66 * dfdy(z, n, i) * dfdy(z, n,
                                                                                                             j)) * np.linalg.det(
                        J(z, n))

                    def EKuv(z, n): return (C12 * dfdx(z, n, i) * dfdy(z, n, j) + C66 * dfdy(z, n, i) * dfdx(z, n,
                                                                                                             j)) * np.linalg.det(
                        J(z, n))

                    def EKvu(z, n): return (C12 * dfdy(z, n, i) * dfdx(z, n, j) + C66 * dfdx(z, n, i) * dfdy(z, n,
                                                                                                             j)) * np.linalg.det(
                        J(z, n))

                    def EKvv(z, n): return (C11 * dfdy(z, n, i) * dfdy(z, n, j) + C66 * dfdx(z, n, i) * dfdx(z, n,
                                                                                                             j)) * np.linalg.det(
                        J(z, n))

                    Kuu[i, j] = e.intGauss2D(gauss, EKuu)
                    Kuv[i, j] = e.intGauss2D(gauss, EKuv)
                    Kvu[i, j] = e.intGauss2D(gauss, EKvu)
                    Kvv[i, j] = e.intGauss2D(gauss, EKvv)

                def EFu(z, n): return (Fx(x(z, n), y(z, n))) * \
                    np.linalg.det(J(z, n))

                def EFv(z, n): return (Fy(x(z, n), y(z, n))) * \
                    np.linalg.det(J(z, n))

                Fu[i] = e.intGauss2D(gauss, EFu)
                Fv[i] = e.intGauss2D(gauss, EFv)

            MD = np.linspace(0, 2 * n - 1, 2 * n).reshape([2, n]).astype(int)
            K[np.ix_(MD[0], MD[0])] = Kuu
            K[np.ix_(MD[0], MD[1])] = Kuv
            K[np.ix_(MD[1], MD[0])] = Kvu
            K[np.ix_(MD[1], MD[1])] = Kvv
            np.savetxt(f'KL_{count}.csv', t*K, delimiter=',')

            F[np.ix_(MD[0])] = Fu
            F[np.ix_(MD[1])] = Fv
            Knl = []
            countnl = 0
            dzpsis = e.dzpsi
            dnpsis = e.dnpsi
            psi = e.psi
            for p in e.elementosnl:
                enl = this.elementos[p]
                countnl += 1
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

                def A(X, XP, l=0.1, t=0.5):
                    l0 = (1)/(2*np.pi*l**2*t)
                    return l0 * np.exp(-distancia(x(X[0], X[1]), y(X[0], X[1]), xnl(XP[0], XP[1]), ynl(XP[0], XP[1])) / l)
                dzpsisnl = enl.dzpsi
                dnpsisnl = enl.dnpsi
                psinl = enl.psi
                for i in range(n):
                    for j in range(nnl):

                        def FUNCIONES(z, n, znl, nnl):

                            jacobianonl = Jnl(znl, nnl)
                            jacobiano_nl = _Jnl(znl, nnl)
                            jacobiano_ = _J(z, n)
                            jacobiano = J(z, n)

                            detjacnl = np.linalg.det(jacobianonl)
                            detjac = np.linalg.det(jacobiano)
                            dz_i = dzpsis[i](z, n)
                            dn_i = dnpsis[i](z, n)
                            psis_i = psi[i](z, n)

                            dfdx_i = dz_i * \
                                jacobiano_[0][0] + dn_i * jacobiano_[0][1]
                            dfdy_i = dz_i * \
                                jacobiano_[1][0] + dn_i * jacobiano_[1][1]

                            dznl_j = dzpsisnl[j](znl, nnl)
                            dnnl_j = dnpsisnl[j](znl, nnl)
                            psisnl_j = psinl[j](znl, nnl)

                            dfdxnl_j = dznl_j * \
                                jacobiano_nl[0][0] + \
                                dnnl_j * jacobiano_nl[0][1]
                            dfdynl_j = dznl_j * \
                                jacobiano_nl[1][0] + \
                                dnnl_j * jacobiano_nl[1][1]

                            AZN = A([z, n], [znl, nnl])

                            EKuunl = AZN * \
                                (C11 * dfdx_i * dfdxnl_j + C66 *
                                 dfdy_i * dfdynl_j) * detjac * detjacnl
                            EKuvnl = AZN * \
                                (C12 * dfdx_i * dfdynl_j + C66 *
                                 dfdy_i * dfdxnl_j) * detjac * detjacnl
                            EKvunl = AZN * \
                                (C12 * dfdy_i * dfdxnl_j + C66 *
                                 dfdx_i * dfdynl_j) * detjac * detjacnl
                            EKvvnl = AZN * \
                                (C11 * dfdy_i * dfdynl_j + C66 *
                                 dfdx_i * dfdxnl_j) * detjac * detjacnl
                            return np.array([EKuunl, EKuvnl, EKvunl, EKvvnl])
                        matrices = e.intGauss4DNL(
                            gaussnl, FUNCIONES, shape=[4, 1])
                        Kuunl[i, j] = matrices[0]
                        Kuvnl[i, j] = matrices[1]
                        Kvunl[i, j] = matrices[2]
                        Kvvnl[i, j] = matrices[3]

                Knlmn[np.ix_(MD[0], MD[0])] = Kuunl
                Knlmn[np.ix_(MD[0], MD[1])] = Kuvnl
                Knlmn[np.ix_(MD[1], MD[0])] = Kvunl
                Knlmn[np.ix_(MD[1], MD[1])] = Kvvnl
                Knl.append(Knlmn)
                progressbar(countnl, len(e.elementosnl),
                            prefix="Integrando elementos No Locales   ", size=50)

            e.determinarMatrices(t*K, F, Q)
            e.KNLS = t**2*np.array(Knl)
            KNLS__ = []

            for i in range(len(e.KNLS)):
                KNLS__.append(e.KNLS[i].T.reshape([1, int(4*n**2)])[0])
            np.savetxt(f'KNLS_{count}.csv', KNLS__, delimiter=',')

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

            e.coordsNuevas = np.array(
                coords) + Ue.T.reshape([2, len(e.psis(0, 0))]).T * mult
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

            progressbar(count, len(this.elementos),
                        prefix="Graficando Solucion ", size=50)
            ax.plot(X_0, Y_0, 'o--', color='gray', zorder=0)
            ax.plot(X_F, Y_F, 'o-', color='black', zorder=3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Solución de elementos finitos')
        ax.legend(['Situacion Inicial', 'Deformada'])

    def defUnitariaX(this, figsize):
        count = 0
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(projection='3d')
        markersize = 2
        cmap = 'magma'
        linewidth = 2
        xtotal = []
        ytotal = []
        ztotal = []
        count = 0
        print('Graficando deformaciones...')
        for e in tqdm(this.elementos, unit='Elemento'):
            count += 1
            x, y, u = postProcesoX(e, this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
        surf = ax.plot_trisurf(xtotal, ytotal, ztotal, cmap=cmap, zorder=1)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel(r'$\varepsilon x$')
        ax.set_title(r'$\frac{\partial U}{\partial X}$')
        print('Deformaciones graficadas')
        # plt.show()
        return xtotal, ytotal, ztotal

    def defUnitariaXMOD(this, x0, xf, xx):
        count = 0
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        markersize = 2
        cmap = 'magma'
        linewidth = 2
        xtotal = []
        ytotal = []
        ztotal = []
        count = 0
        print('Graficando deformaciones...')
        for e in tqdm(this.elementos, unit='Elemento'):
            count += 1
            x, y, u = postProcesoX(e, this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
        surf = ax.plot_trisurf(xtotal, ytotal, ztotal,
                               cmap=cmap, zorder=1, alpha=0.3)
        _Y = np.linspace(x0, xf, 100).tolist()
        Z = []
        Y = []
        Zlocal = []
        print('Graficando perfil')
        for y in tqdm(_Y):
            for i in range(len(this.elementos)):
                e = this.elementos[i]
                if e.estaDentro(y, xx):
                    U = postProcesoXNodos(e, this.U)
                    Y.append(y)
                    Z.append(e.darSolucionXY(U, y, xx, n=100)[0])
                    break
        ax.plot(Y, [xx]*len(Y), Z)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        print('Deformaciones graficadas')

    def defUnitariaY(this, figsize):
        count = 0
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
            count += 1
            x, y, u = postProcesoY(e, this.U)
            xtotal.extend(x)
            ytotal.extend(y)
            ztotal.extend(u)
        surf = ax.plot_trisurf(xtotal, ytotal, ztotal, cmap=cmap, zorder=1)
        cbar = fig.colorbar(surf)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel(r'$\varepsilon x$')
        ax.set_title(r'$\frac{\partial U}{\partial Y}$')

    def darDeformacionXY(this, U, x, y, n=100):
        this.Ue = U[np.ix_(this.gdl)]
        this.U = lambda z, n: this.Ue.T[0] @ this.psis(z, n)
        zeta, eta = this.mappingInverso(x, y, n).reshape(1, 2)[0].tolist()
        return this.U(zeta, eta)

    def ensamblar(this):
        z1 = this.z1
        z2 = 1 - z1
        print('Ensamblando...')
        for e in tqdm(this.elementos, unit='Elemento'):
            this.K[np.ix_(e.gdl, e.gdl)] += e.Ke*z1
            for i, enl in enumerate(np.array(this.elementos)[e.elementosnl]):
                this.K[np.ix_(e.gdl, enl.gdl)] += e.KNLS[i]*z2
            this.F[np.ix_(e.gdl)] += e.Fe
            this.Q[np.ix_(e.gdl)] += e.Qe
        this._K = np.copy(this.K)
        print('Sistema de ecuaciones ensamblado')

    def perfilX(Objeto_FEM, x, yi=0.00016, yf=0.0004, x0=0, xf=1, filename='', acum=False, label=''):
        _Y = np.linspace(x0, xf, 100).tolist()
        Z = []
        Y = []
        Zlocal = []
        print('Graficando perfil')
        for y in tqdm(_Y):
            for i in range(len(Objeto_FEM.elementos)):
                e = Objeto_FEM.elementos[i]
                if e.estaDentro(x, y):
                    print(i, x, y)
                    U = postProcesoXNodos(e, Objeto_FEM.U)
                    try:
                        Z.append(e.darSolucionXY(U, x, y, n=100))
                        Y.append(y)
                    except:
                        print('No se pudo graficár la coordenada',
                              x, y, "del elemento", i)
                    break
        plt.plot(Y, Z, label=label)
        print('Perfil graficado')

        if acum:
            plt.grid()
            plt.legend()
            plt.xlabel('y')
            if not filename == '':
                plt.savefig(filename, transparent=True)
            # plt.show()

    def perfilY(Objeto_FEM, y, yi=0.00016, yf=0.0004, x0=0, xf=50, filename='', acum=False, label=''):
        _X = np.linspace(x0, xf, 100).tolist()
        Z = []
        X = []
        Zlocal = []
        print('Graficando perfil')
        for x in tqdm(_X):
            for i in range(len(Objeto_FEM.elementos)):
                e = Objeto_FEM.elementos[i]
                if e.estaDentro(x, y):
                    print(i, x, y)
                    U = postProcesoXNodos(e, Objeto_FEM.U)
                    try:
                        Z.append(e.darSolucionXY(U, x, y, n=100))
                        X.append(x)
                    except:
                        print('No se pudo graficár la coordenada',
                              x, y, "del elemento", i)
                    break
        plt.plot(X, Z, label=label)
        print('Perfil graficado')

        if acum:
            plt.legend()
            plt.grid()
            plt.xlabel('x')
            if not filename == '':
                plt.savefig(filename, transparent=True)
            # plt.show()
        return X, Z

    def moduloIntegrador(this, enmallado, npg, E, V, t, l, RUTA_M, tfa=1):
        subprocess.run("index.exe"+" " + enmallado + " "+format(int(npg))+" "+format(E)+" "+format(
            V)+" "+format(t)+" "+format(l)+" "+RUTA_M+" "+format(int(tfa)), shell=True, check=True)

    def moduloIntegradorParalelo(this, enmallado, npg, E, V, t, l, RUTA_M, tfa=1):
        subprocess.Popen("index.exe"+" " + enmallado + " "+format(int(npg))+" "+format(
            E)+" "+format(V)+" "+format(t)+" "+format(l)+" "+RUTA_M+" "+format(int(tfa)))


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
    return np.array(u).reshape([len(Z), 1])


def postProcesoY(this, U):
    this.Ue = U[np.ix_(this.gdl)]
    this._Ue = this.Ue.T[0].tolist()
    this._Ue.append(this.Ue[0][0])
    Z = this._dominioNaturalZ
    N = this._dominioNaturalN
    x = []
    y = []
    u = []
    for z, n in zip(Z, N):
        x.append(this.Tx(z, n)[0])
        y.append(this.Ty(z, n)[0])
        u.append(grad(this, z, n)[1][0])
    return x, y, u


def grad(this, z, n):
    dz = this.dzpsis(z, n)
    dn = this.dnpsis(z, n)
    result = []
    _j = this._J(z, n)
    for i in range(len(dz)):
        result.append(_j @ np.array([[dz[i][0]], [dn[i][0]]]))
    result = np.array(result)
    n = len(result)
    U = np.linspace(0, n-1, n).astype(int)
    V = U*2+1
    dx = (this.Ue[np.ix_(U)].T @ result[:, 0])[0]
    dy = (this.Ue[np.ix_(U)].T @ result[:, 1])[0]
    return np.array([dx, dy])
