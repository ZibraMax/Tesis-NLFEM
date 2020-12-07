import numpy as np
from ..Utils import isBetween
import matplotlib.pyplot as plt

class Geometria:
    def __init__(this, vertices, diccionarios, gdls, tipos, segmentos=[]):
        this.Areas = []
        this.vertices = vertices
        this.diccionarios = diccionarios
        this.gdls = gdls
        this.tipos = tipos
        this.segmentos = segmentos
        this._diccionarios = []
        this.cbe = []
        this.cbn = []
        this.centroideYArea()

    def centroideYArea(this):
        for i,e in enumerate(this.diccionarios):
            if this.tipos[i]=='T2V':
                f4 = e[5]
                f5 = e[3]
                f6 = e[4]
                e[3] = f4
                e[4] = f5
                e[5] = f6
                this._diccionarios.append([e[0],e[1],e[2]])
            elif this.tipos[i]=='C1V':
                this._diccionarios.append([e[0],e[1],e[2],e[3]])
        this.Centroides = []
        for i, e in enumerate(this._diccionarios):
            coords = np.array(this.gdls)[np.ix_(e)]
            coords = np.array(coords.tolist() + [coords[0].tolist()])
            area = 0
            cx = 0
            cy = 0
            for j in range(len(coords)-1):
                area += coords[j][0]*coords[j+1][1]-coords[j+1][0]*coords[j][1]
                mult = (coords[j][0]*coords[j+1][1]-coords[j+1][0]*coords[j][1])
                cx += (coords[j][0]+coords[j+1][0])*mult
                cy += (coords[j][1]+coords[j+1][1])*mult
            this.Areas.append(np.abs(area/2))
            this.Centroides.append([cx/3/area,cy/3/area])
    def generarSegmentosDesdeCoordenadas(this,p0,p1):
        masCercano1 = None
        d1 = np.Inf
        masCercano2 = None
        d2 = np.Inf
        for i,gdl in enumerate(this.gdls):
            r1 = np.sqrt((p0[0]-gdl[0])**2+(p0[1]-gdl[1])**2)
            r2 = np.sqrt((p1[0]-gdl[0])**2+(p1[1]-gdl[1])**2)
            if r1 < d1:
                d1 = r1
                masCercano1 = i
            if r2 < d2:
                d2 = r2
                masCercano2 = i
        this.segmentos.append([masCercano1,masCercano2])
    def generarCBdesdeCoordenada(this,x,y,valor=0):
        masCercano1 = None
        d1 = np.Inf
        for i,gdl in enumerate(this.gdls):
            r1 = np.sqrt((x-gdl[0])**2+(y-gdl[1])**2)
            if r1 < d1:
                d1 = r1
                masCercano1 = i
        return [[i,valor]]
    def darNodosCB(this, segmento):
        a = []
        ps = np.array(this.gdls)[this.segmentos[segmento]].tolist()
        for i, p in enumerate(this.gdls):
            if isBetween(ps[0], ps[1], p):
                a.append(i)
        return np.array(a)

    def generarCBdesdeBorde(this, borde, valor=0):
        cb = []
        nodos = this.darNodosCB(borde)
        cbe = np.zeros([len(nodos), 2])
        cbe[:, 0] = nodos
        cbe[:, 1] = valor
        cb += cbe.tolist()
        return cb
    def definirTodasCondiciones(this):
        for s in range(len(this.segmentos)):
            this.cbe += this.generarCBdesdeBorde(s)
    def generarCBdesdeBordeX(this, borde, valor=0):
        cb = []
        nodos = this.darNodosCB(borde)
        cbe = np.zeros([len(nodos), 2])
        cbe[:, 0] = nodos * 2
        cbe[:, 1] = valor
        cb += cbe.tolist()
        return cb

    def generarCBdesdeBordeY(this, borde, valor=0):
        cb = []
        nodos = this.darNodosCB(borde)
        cbe = np.zeros([len(nodos), 2])
        cbe[:, 0] = nodos * 2 + 1
        cbe[:, 1] = valor
        cb += cbe.tolist()
        return cb
    def dibujarse(this,texto=10,bolita=0,figsize=[17,10]):

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()

        ax.axes.set_aspect('equal')

        for i, e in enumerate(this._diccionarios):
            coords = np.array(this.gdls)[np.ix_(e)]
            coords = np.array(coords.tolist() + [coords[0].tolist()])
            X = coords[:, 0]
            Y = coords[:, 1]
            ax.plot(X, Y, 'o-', color='black', zorder=-10)
            cx = this.Centroides[i][0]
            cy = this.Centroides[i][1]
            ax.plot(cx, cy, 'o', markersize=texto + bolita, color='yellow')
            ax.annotate(format(i), [cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        try:
            verts = this.gdls
            segs = this.segmentos
            for i,seg in enumerate(segs):
                x0, y0 = verts[int(seg[0])]
                x1, y1 = verts[int(seg[1])]

                ax.fill(
                    [x0, x1],
                    [y0, y1],
                    facecolor='none',
                    edgecolor='b',
                    linewidth=3,
                    zorder=0,
                )
                cx = (x0+x1)*0.5
                cy = (y0+y1)*0.5
                ax.plot(cx, cy, 'o', markersize=texto + bolita, color='pink')
                ax.annotate(format(i), [cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        except:
            pass
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Dominio')

        gdls = np.array(this.gdls)

        labels = np.linspace(0, gdls.shape[0] - 1, gdls.shape[0]).astype(int)

        ax.plot(gdls[:, 0], gdls[:, 1], 'o', markersize=texto+bolita, color='gray')

        for p, l in zip(gdls, labels):
            ax.annotate(l, p, size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
