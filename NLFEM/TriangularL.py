from .integradorTriangular import *
class TriangularL(integradorTriangular):
    def __init__(this,coords,gdl=None,gauss=4):
        if len(coords) == 3:
            this.ZNatural = [0,1,0]
            this.NNatural = [0,0,1]
            this.psis = lambda z,n: np.array([[1-z-n],[z],[n]])
            this.dzpsis = lambda z,n: np.array([[-1],[1],[0]])
            this.dnpsis = lambda z,n: np.array([[-1],[0],[1]])
        else:
            raise Exception('Error: Se esta intentando crear un elemento triangualar que no tiene 3 coordenadas. '+'Recuerde que se necesitan de 3 corrdenadas en orden contrario de las manecillas del reloj '+'sin repetir la primera coordenada')
        super().__init__(coords=coords,gdl=gdl,gauss=gauss)
        this._coords = coords