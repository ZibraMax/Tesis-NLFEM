#include <vector>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Elemento.h"
using namespace std;

class TriangularC : public Elemento {
	public:
		TriangularC(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2],vector<double> ppgauss,vector<double> pwgauss,double pt):Elemento(pgdls,pnolocales,COORDENADAS,ppgauss,pwgauss,pt) {
			psi = {
			[](double z, double n){return 2.0*(z+n-1.0)*(z+n-1.0/2.0);},
			[](double z, double n){return 2.0*z*(z-1.0/2.0);},
			[](double z, double n){return 2.0*n*(n-1.0/2.0);},
			[](double z, double n){return -4.0*(z+n-1.0)*(z);},
			[](double z, double n){return 4.0*z*n;},
			[](double z, double n){return -4.0*n*(z+n-1.0);}};

			dzpsi = {
			[](double z, double n){return 4.0*z+4.0*n-3.0;},
			[](double z, double n){return 4.0*z-1.0;},
			[](double z, double n){return 0.0;},
			[](double z, double n){return -8.0*z-4.0*(n-1.0);},
			[](double z, double n){return 4.0*n;},
			[](double z, double n){return -4.0*n;}};

			dnpsi = {
			[](double z, double n){return 4.0*n+4.0*z-3.0;},
			[](double z, double n){return 0.0;},
			[](double z, double n){return 4.0*n-1.0;},
			[](double z, double n){return -4.0*z;},
			[](double z, double n){return 4.0*z;},
			[](double z, double n){return -8.0*n-4.0*z+4.0;}};
			double A0 = 1.0/3.0;
	        double A1 = 0.059715871789770;
	        double A2 = 0.797426985353087;
	        double B1 = 0.470142064105115;
	        double B2 = 0.101286507323456;
	        double W0 = 0.1125;
	        double W1 = 0.066197076394253;
	        double W2 = 0.062969590272413;
	        vector<double> X = {A0,A1,B1,B1,B2,B2,A2};
	        vector<double> Y = {A0,B1,A1,B1,A2,B2,B2};
	        PESOS = {W0,W1,W1,W1,W2,W2,W2};
			for (int i = 0; i < X.size(); ++i) {
				vector<double> puntoi = {X[i], Y[i]};
				PUNTOS.push_back(puntoi);
			}
		}
};