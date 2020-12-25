#include <vector>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Elemento.h"
using namespace std;

class Serendipity : public Elemento {
	public:
		Serendipity(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2],vector<double> ppgauss,vector<double> pwgauss,double pt):Elemento(pgdls,pnolocales,COORDENADAS,ppgauss,pwgauss,pt) {
			psi = {
			[](double z, double n){return 1.0/4.0*(1.0-z)*(1.0-n)*(-1.0-z-n);},
			[](double z, double n){return 1.0/4.0*(1.0+z)*(1.0-n)*(-1.0+z-n);},
			[](double z, double n){return 1.0/4.0*(1.0+z)*(1.0+n)*(-1.0+z+n);},
			[](double z, double n){return 1.0/4.0*(1.0-z)*(1.0+n)*(-1.0-z+n);},
			[](double z, double n){return 1.0/2.0*(1.0-z*z)*(1.0-n);},
			[](double z, double n){return 1.0/2.0*(1.0+z)*(1.0-n*n);},
			[](double z, double n){return 1.0/2.0*(1.0-z*z)*(1.0+n);},
			[](double z, double n){return 1.0/2.0*(1.0-z)*(1.0-n*n);}};

			dzpsi = {
			[](double z, double n){return -1.0/4.0*(n-1.0)*(2.0*z+n);},
			[](double z, double n){return -1.0/4.0*(n-1.0)*(2.0*z-n);},
			[](double z, double n){return 1.0/4.0*(n+1.0)*(2.0*z+n);},
			[](double z, double n){return 1.0/4.0*(n+1.0)*(2.0*z-n);},
			[](double z, double n){return (n-1.0)*z;},
			[](double z, double n){return -1.0/2.0*(n*n-1.0);},
			[](double z, double n){return -(n+1.0)*z;},
			[](double z, double n){return 1.0/2.0*(n*n-1.0);}};

			dnpsi = {
			[](double z, double n){return -1.0/4.0*(z-1.0)*(2.0*n+z);},
			[](double z, double n){return 1.0/4.0*(z+1.0)*(2.0*n-z);},
			[](double z, double n){return 1.0/4.0*(z+1.0)*(2.0*n+z);},
			[](double z, double n){return -1.0/4.0*(z-1.0)*(2.0*n-z);},
			[](double z, double n){return 1.0/2.0*(z*z-1.0);},
			[](double z, double n){return -n*(z+1.0);},
			[](double z, double n){return -1.0/2.0*(z*z-1.0);},
			[](double z, double n){return n*(z-1.0);}};

			for (int i = 0; i < ppgauss.size(); ++i) {
				for (int j = 0; j < ppgauss.size(); ++j) {
					vector<double> puntoi = {ppgauss[i], ppgauss[j]};
					PUNTOS.push_back(puntoi);
					PESOS.push_back(pwgauss[i]*pwgauss[j]);
				}
			}
		}
		vector<double> psis(double z, double n) {
			vector<double> foo{
				1.0/4.0*(1.0-z)*(1.0-n)*(-1.0-z-n),
				1.0/4.0*(1.0+z)*(1.0-n)*(-1.0+z-n),
				1.0/4.0*(1.0+z)*(1.0+n)*(-1.0+z+n),
				1.0/4.0*(1.0-z)*(1.0+n)*(-1.0-z+n),
				1.0/2.0*(1.0-z*z)*(1.0-n),
				1.0/2.0*(1.0+z)*(1.0-n*n),
				1.0/2.0*(1.0-z*z)*(1.0+n),
				1.0/2.0*(1.0-z)*(1.0-n*n)
				};
			return foo;
		}
		vector<double> dzpsis(double z, double n) {
			vector<double> foo{
				-1.0/4.0*(n-1.0)*(2.0*z+n),
				-1.0/4.0*(n-1.0)*(2.0*z-n),
				1.0/4.0*(n+1.0)*(2.0*z+n),
				1.0/4.0*(n+1.0)*(2.0*z-n),
				(n-1.0)*z,
				-1.0/2.0*(n*n-1.0),
				-(n+1.0)*z,
				1.0/2.0*(n*n-1.0)
				};
			return foo;
		}
		vector<double> dnpsis(double z, double n) {
			vector<double> foo{
				-1.0/4.0*(z-1.0)*(2.0*n+z),
				1.0/4.0*(z+1.0)*(2.0*n-z),
				1.0/4.0*(z+1.0)*(2.0*n+z),
				-1.0/4.0*(z-1.0)*(2.0*n-z),
				1.0/2.0*(z*z-1.0),
				-n*(z+1.0),
				-1.0/2.0*(z*z-1.0),
				n*(z-1.0)
				};
			return foo;
		}
};