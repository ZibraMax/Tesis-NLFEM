from .integradorRectangular import *
class SerendipityC(integradorRectangular):
    def __init__(this,coords,gdl=None,gauss=4):
        if len(coords) == 8:
            this.psis = lambda z,n: np.array([[1/4*(1-z)*(1-n)*(-1-z-n)],[1/4*(1+z)*(1-n)*(-1+z-n)],[1/4*(1+z)*(1+n)*(-1+z+n)],[1/4*(1-z)*(1+n)*(-1-z+n)],[1/2*(1-z**2)*(1-n)],[1/2*(1+z)*(1-n**2)],[1/2*(1-z**2)*(1+n)],[1/2*(1-z)*(1-n**2)]])
            this.psi = np.array([lambda z,n: 1/4*(1-z)*(1-n)*(-1-z-n),
            	lambda z,n: 1/4*(1+z)*(1-n)*(-1+z-n),
            	lambda z,n: 1/4*(1+z)*(1+n)*(-1+z+n),
            	lambda z,n: 1/4*(1-z)*(1+n)*(-1-z+n),
            	lambda z,n: 1/2*(1-z**2)*(1-n),
            	lambda z,n: 1/2*(1+z)*(1-n**2),
            	lambda z,n: 1/2*(1-z**2)*(1+n),
            	lambda z,n: 1/2*(1-z)*(1-n**2)])

            this.dzpsis = lambda z,n: np.array([[-1/4*(n-1)*(2*z+n)],[-1/4*(n-1)*(2*z-n)],[1/4*(n+1)*(2*z+n)],[1/4*(n+1)*(2*z-n)],[(n-1)*z],[-1/2*(n**2-1)],[-(n+1)*z],[1/2*(n**2-1)]])
            this.dzpsi = np.array([lambda z,n: -1/4*(n-1)*(2*z+n),
            	lambda z,n: -1/4*(n-1)*(2*z-n),
            	lambda z,n: 1/4*(n+1)*(2*z+n),
            	lambda z,n: 1/4*(n+1)*(2*z-n),
            	lambda z,n: (n-1)*z,
            	lambda z,n: -1/2*(n**2-1),
            	lambda z,n: -(n+1)*z,
            	lambda z,n: 1/2*(n**2-1)])

            this.dnpsis = lambda z,n: np.array([[-1/4*(z-1)*(2*n+z)],[1/4*(z+1)*(2*n-z)],[1/4*(z+1)*(2*n+z)],[-1/4*(z-1)*(2*n-z)],[1/2*(z**2-1)],[-n*(z+1)],[-1/2*(z**2-1)],[n*(z-1)]])
            this.dnpsi = np.array([lambda z,n: -1/4*(z-1)*(2*n+z),
            	lambda z,n: 1/4*(z+1)*(2*n-z),
            	lambda z,n: 1/4*(z+1)*(2*n+z),
            	lambda z,n: -1/4*(z-1)*(2*n-z),
            	lambda z,n: 1/2*(z**2-1),
            	lambda z,n: -n*(z+1),
            	lambda z,n: -1/2*(z**2-1),
            	lambda z,n: n*(z-1)])

        else:
            raise Exception('Error: Se esta intentando crear un elemento que no tiene 8 coordenadas. '+'Recuerde que se necesitan de 8 corrdenadas en sentido contrario de las manecillas del reloj '+'sin repetir la primera coordenada')
        super().__init__(coords=coords,gdl=gdl,gauss=gauss)