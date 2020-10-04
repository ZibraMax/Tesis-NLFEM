import triangle as tr

import numpy as np
import matplotlib.pyplot as plt

from .Geometria import *
import pandas as pd

class Delaunay1V(Geometria):
    def __init__(this, vertices, params, plot=False,texto=10,bolita=0):
        this.params = params
        this.cbe = []
        this.cbn = []
        this.seg = []
        for i in range(len(vertices)-1):
            this.seg.append([i,i+1])
        this.seg.append([i+1,0])
        this.original = dict(vertices=np.array(vertices),segments=np.array(this.seg))
        this.triangular = tr.triangulate(this.original,this.params)
        diccionarios = this.triangular['triangles'].tolist()
        tipos = np.zeros([len(diccionarios)]).astype(str)
        if 'o2' in params:
            tipos[:] = 'T2V'
        else:
            tipos[:] = 'T1V'
        gdls = this.triangular['vertices'].tolist()
        super().__init__(vertices,diccionarios,gdls,tipos,segmentos=this.seg)
        if plot:
            this.dibujarse(texto)

def _strdelaunay(constrained=True,delaunay=True,a=None,q=None,o=1):
    p = ''
    if o==2:
        o = '-o2'
    else:
        o = ''
    if constrained:
        p = 'p'
    if a == None:
        a = ''
    else:
        a = 'a'+format(a)
    D = ''
    if delaunay:
        D = 'D'
    if q == None:
        q=''
    else:
        if type(q) == int:
            if q > 35:
                raise "No se puede crear una triangulacion con angulos menores a 35 grados"
        q = 'q'+format(q)
    return p+a+D+q+'i'+o
def generarGeometriaDesdeTriangulacion(triang):
    count = 0
    for i in triang['segments']:
        if count > 1:
            if np.sum(np.isin(np.array(cbe)[:,0], i[0]))<1:
                cbe.append([i[0],0])
            if np.sum(np.isin(np.array(cbe)[:,0], i[1]))<1:
                cbe.append([i[1],0])
        else:
            cbe.append([i[0],0])
        if np.sum(np.isin(np.array(cbe)[:,0], i[1]))<1:
                cbe.append([i[1],0])
        count+=1
    dic = triang['triangles'].tolist()
    tipos = np.zeros([len(dic)]).astype(str)
    tipos[:] = 'T1V'
    g = Geometria([-1],dic,triang['vertices'].tolist(),tipos,segmentos=triang['segments'])
    g.triangular = triang
    return g
    