import numpy as np
from IPython.display import clear_output

class FEM:
    def __init__(this,ngdl):
        this.elementos = []
        this.n = ngdl
        this.K = np.zeros([this.n,this.n])
        this.F = np.zeros([this.n,1])
        this.Q = np.zeros([this.n,1])
        this.U = np.zeros([this.n,1])
        this.S = np.zeros([this.n,1])
        this._K = np.copy(this.K)
        this.cbe = []
        this.cbn = []
        this.geometria = None
    def importarGeometria(this,geometria):
        this.geometria = geometria
    def ensamblar(this):
        for e in this.elementos:
            this.K[np.ix_(e.gdl,e.gdl)] += e.Ke
            this.F[np.ix_(e.gdl)] += e.Fe
            this.Q[np.ix_(e.gdl)] += e.Qe
        this._K = np.copy(this.K)
    def condicionesFrontera(this):
        for i in this.cbn:
            this.Q[int(i[0])] = i[1]
        for i in this.cbe:
            ui = np.zeros([this.n, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(this.K, ui)
            this.S = this.S - vv
            this.K[int(i[0]), :] = 0
            this.K[:, int(i[0])] = 0
            this.K[int(i[0]), int(i[0])] = 1
        this.S = this.S + this.F + this.Q
        for i in this.cbe:
            this.S[int(i[0])] = i[1]
    def solucionarSistemaEcuaciones(this):
        this.U = np.linalg.solve(this.K,this.S)
        for e in this.elementos:
            e.Ue = this.U[np.ix_(e.gdl)]
    def solucionar(this,plot=True,figsize=[14,12],cmap='magma',markersize=2,linewidth=2,mask=None, **kargs):
        this.generarElementos()
        this.calcularMatrices(**kargs)
        print('Ensamblando sistema de ecuaciones')
        this.ensamblar()
        print('Definiendo condiciones deborde')
        this.condicionesFrontera(this.cbe,this.cbn)
        print('Solucionando sistema de ecuaciones')
        this.solucionarSistemaEcuaciones()
        if plot:
            print('Graficando...')
            clear_output(wait=True)
            this.graficarSolucion(figsize,cmap,mask=mask)
    def definirCondicionesDeBorde(this,cbe=[],cbn=[]):
        this.cbe=cbe
        this.cbn=cbn
def progressbar(i, length, prefix="", size=60):
    def show(j):
        clear_output(wait=True)
        x = int(size*j/length)
        print("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, length))
    show(i)