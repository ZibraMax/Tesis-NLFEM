#include <vector>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Elemento.h"
using namespace std;

class Rectangular : public Elemento {
	public:
		Rectangular(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2],vector<double> ppgauss,vector<double> pwgauss,double pt):Elemento(pgdls,pnolocales,COORDENADAS,ppgauss,pwgauss,pt) {
			psi = {
			[](double z, double n){return 1.0/4.0*(1.0-z)*(1.0-n);},
			[](double z, double n){return 1.0/4.0*(1.0+z)*(1.0-n);},
			[](double z, double n){return 1.0/4.0*(1.0+z)*(1.0+n);},
			[](double z, double n){return 1.0/4.0*(1.0-z)*(1.0+n);}};

			dzpsi = {
			[](double z, double n){return 1.0/4.0*(n-1.0);},
			[](double z, double n){return -1.0/4.0*(n-1.0);},
			[](double z, double n){return 1.0/4.0*(n+1.0);},
			[](double z, double n){return -1.0/4.0*(1.0+n);}};

			dnpsi = {
			[](double z, double n){return 1.0/4.0*(z-1.0);},
			[](double z, double n){return -1.0/4.0*(z+1.0);},
			[](double z, double n){return 1.0/4.0*(1.0+z);},
			[](double z, double n){return 1.0/4.0*(1.0-z);}};

			for (int i = 0; i < ppgauss.size(); ++i) {
				for (int j = 0; j < ppgauss.size(); ++j) {
					vector<double> puntoi = {ppgauss[i], ppgauss[j]};
					PUNTOS.push_back(puntoi);
					PESOS.push_back(pwgauss[i]*pwgauss[j]);
				}
			}
		}
};