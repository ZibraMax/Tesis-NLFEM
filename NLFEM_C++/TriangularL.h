#include <vector>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Elemento.h"
using namespace std;

class TriangularL : public Elemento {
	public:
		TriangularL(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2],vector<double> ppgauss,vector<double> pwgauss,double pt):Elemento(pgdls,pnolocales,COORDENADAS,ppgauss,pwgauss,pt) {
			psi = {
			[](double z, double n){return 1.0-z-n;},
			[](double z, double n){return z;},
			[](double z, double n){return n;}};

			dzpsi = {
			[](double z, double n){return -1;},
			[](double z, double n){return 1;},
			[](double z, double n){return 0;}};

			dnpsi = {
			[](double z, double n){return -1;},
			[](double z, double n){return 0;},
			[](double z, double n){return 1;}};

			for (int i = 0; i < ppgauss.size(); ++i) {
				for (int j = 0; j < ppgauss.size(); ++j) {
					vector<double> puntoi = {ppgauss[i], ppgauss[j]};
					PUNTOS.push_back(puntoi);
					PESOS.push_back(pwgauss[i]*pwgauss[j]);
				}
			}
		}
};