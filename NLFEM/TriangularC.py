from .integradorTriangular import *
class TriangularC(integradorTriangular):
    def __init__(this,coords,gdl=None,gauss=4):
        if len(coords) == 6:
            this.ZNatural = [0,1,0,0.5,0.5,0]
            this.NNatural = [0,0,1,0,0.5,0.5]
            this.psi = np.array([lambda z,n: 2*(z+n-1)*(z+n-1/2),
                                 lambda z,n: 2*z*(z-1/2),
                                 lambda z,n: 2*n*(n-1/2),
                                 lambda z,n: -4*(z+n-1)*(z),
                                 lambda z,n: 4*z*n,
                                 lambda z,n: -4*n*(z+n-1)])
            this.psis = lambda z,n: np.array([[2*(z+n-1)*(z+n-1/2)],
                                 [2*z*(z-1/2)],
                                 [2*n*(n-1/2)],
                                 [-4*(z+n-1)*(z)],
                                 [4*z*n],
                                 [-4*n*(z+n-1)]])
            this.dzpsi = np.array([lambda z,n: 4*z+4*n-3,
                                   lambda z,n: 4*z-1,
                                   lambda z,n: 0,
                                   lambda z,n: -8*z-4*(n-1),
                                   lambda z,n: 4*n,
                                   lambda z,n: -4*n])

            this.dzpsis =  lambda z,n: np.array([[4*z+4*n-3],
                                   [4*z-1],
                                   [0],
                                   [-8*z-4*(n-1)],
                                   [4*n],
                                   [-4*n]])

            this.dnpsi = np.array([lambda z,n: 4*n+4*z-3,
                                   lambda z,n: 0,
                                   lambda z,n: 4*n-1,
                                   lambda z,n: -4*z,
                                   lambda z,n: 4*z,
                                   lambda z,n: -8*n-4*z+4])
            this.dnpsis = lambda z,n: np.array([[4*n+4*z-3],
                                   [0],
                                   [4*n-1],
                                   [-4*z],
                                   [4*z],
                                   [-8*n-4*z+4]])
        else:
            raise Exception('Error: Se esta intentando crear un elemento triangular que no tiene 6 coordenadas. '+'Recuerde que se necesitan de 6 coordenadas en orden contrario de las manecillas del reloj '+'sin repetir la primera coordenada')
        super().__init__(coords=coords,gdl=gdl,gauss=gauss)
        this._coords = [coords[0],coords[1],coords[2]]
