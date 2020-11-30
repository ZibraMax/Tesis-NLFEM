#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Eigen/Dense"
#include <chrono> 
using namespace std::chrono;
using namespace Eigen;
using namespace std;

void guardarMatriz(string name, MatrixXd matrix) {
	cout.precision(15);
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
	ofstream Archivo(name.c_str());
    Archivo << matrix.format(CSVFormat);
	Archivo.close();
}

class Serendipity {
	public:
		vector<vector<double>> coords;
		vector<int> gdl;
		vector<int> nolocales;
		vector<double> X;
		vector<double> Y;
		vector<vector<double>> J;
		vector<vector<double>> _J;
		vector<double> PUNTOS;
		vector<double> PESOS;

		double t;
		vector<function<double(double,double)>> psi{
		[](double z, double n){return 1.0/4.0*(1.0-z)*(1.0-n)*(-1.0-z-n);},
		[](double z, double n){return 1.0/4.0*(1.0+z)*(1.0-n)*(-1.0+z-n);},
		[](double z, double n){return 1.0/4.0*(1.0+z)*(1.0+n)*(-1.0+z+n);},
		[](double z, double n){return 1.0/4.0*(1.0-z)*(1.0+n)*(-1.0-z+n);},
		[](double z, double n){return 1.0/2.0*(1.0-z*z)*(1.0-n);},
		[](double z, double n){return 1.0/2.0*(1.0+z)*(1.0-n*n);},
		[](double z, double n){return 1.0/2.0*(1.0-z*z)*(1.0+n);},
		[](double z, double n){return 1.0/2.0*(1.0-z)*(1.0-n*n);}};

		vector<function<double(double,double)>> dzpsi{
		[](double z, double n){return -1.0/4.0*(n-1.0)*(2.0*z+n);},
		[](double z, double n){return -1.0/4.0*(n-1.0)*(2.0*z-n);},
		[](double z, double n){return 1.0/4.0*(n+1.0)*(2.0*z+n);},
		[](double z, double n){return 1.0/4.0*(n+1.0)*(2.0*z-n);},
		[](double z, double n){return (n-1.0)*z;},
		[](double z, double n){return -1.0/2.0*(n*n-1.0);},
		[](double z, double n){return -(n+1.0)*z;},
		[](double z, double n){return 1.0/2.0*(n*n-1.0);}};

		vector<function<double(double,double)>> dnpsi{
		[](double z, double n){return -1.0/4.0*(z-1.0)*(2.0*n+z);},
		[](double z, double n){return 1.0/4.0*(z+1.0)*(2.0*n-z);},
		[](double z, double n){return 1.0/4.0*(z+1.0)*(2.0*n+z);},
		[](double z, double n){return -1.0/4.0*(z-1.0)*(2.0*n-z);},
		[](double z, double n){return 1.0/2.0*(z*z-1.0);},
		[](double z, double n){return -n*(z+1.0);},
		[](double z, double n){return -1.0/2.0*(z*z-1.0);},
		[](double z, double n){return n*(z-1.0);}};

		Serendipity(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2],vector<double> ppgauss,vector<double> pwgauss,double pt) {
			PUNTOS = ppgauss;
			PESOS =pwgauss;
			gdl = pgdls;
			t = pt;
			nolocales = pnolocales;
			coords = darCoordenadas(COORDENADAS,gdl);
			for (int i = 0; i < coords.size(); ++i) {
				X.push_back((double) coords[i][0]);
				Y.push_back((double) coords[i][1]);
			}
		}
		double TX(double z,double n) {
			double res = 0;
			vector<double> psiszn = psis(z,n);
			for (int i = 0; i < psiszn.size(); ++i) {
				res+=X[i]*psiszn[i];
			}
			return res;
		}
		double TY(double z,double n) { 
			double res = 0;
			vector<double> psiszn = psis(z,n);
			for (int i = 0; i < psiszn.size(); ++i) {
				res+=Y[i]*psiszn[i];
			}
			return res;
		}
		vector<vector<double>> darCoordenadas(double COORDENADAS[][2], vector<int> gdls) {
			vector<vector<double>> v;
			for (int i = 0; i < gdls.size(); ++i) {
				vector<double> v1{COORDENADAS[gdls[i]][0],COORDENADAS[gdls[i]][1]};
				v.push_back(v1);
			}
			return v;
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
		vector<double> transfCoordenadas(double z,double n) {
			double dxdz = 0.0;
			double dydz = 0.0;
			double dxdn = 0.0;
			double dydn = 0.0;
			vector<double> dzpsis_zn = dzpsis(z,n);
			vector<double> dnpsis_zn = dnpsis(z,n);
			for (int i = 0; i < coords.size(); ++i) {
				dxdz += X[i]*dzpsis_zn[i];
				dydz += Y[i]*dzpsis_zn[i];
				dxdn += X[i]*dnpsis_zn[i];
				dydn += Y[i]*dnpsis_zn[i];
			}
			J = {{dxdz,dydz},{dxdn,dydn}};
			double determinanteJ = dxdz*dydn-dydz*dxdn;
			_J = {{dydn*1/determinanteJ,-dydz*1/determinanteJ},{-dxdn*1/determinanteJ,dxdz*1/determinanteJ}};
			double determinante_J = _J[0][0]*_J[1][1]-_J[0][1]*_J[1][0];
			vector<double> retorno{determinanteJ,determinante_J};
			return retorno;
		}
		void imprimirse() {
			for (vector<double> vertice: coords) {
				for (double coordenada: vertice) {
					cout<< coordenada << "\t";
				}
				cout << endl;
			}
		}
		void matrizLocal(double C11,double C12,double C66,int nelemento,string rutaBase) {
			int N = gdl.size();
			MatrixXd KUU(N,N);
			MatrixXd KUV(N,N);
			MatrixXd KVU(N,N);
			MatrixXd KVV(N,N);
			for (int gdli = 0; gdli < N; ++gdli) {
				for (int gdlj = 0; gdlj < N; ++gdlj) {
					double KKUU = 0.0;
					double KKUV = 0.0;
					double KKVU = 0.0;
					double KKVV = 0.0;

					for (int i = 0; i < PUNTOS.size(); ++i) {
						double z = PUNTOS[i];
						for (int j = 0; j < PUNTOS.size(); ++j) {
							double n = PUNTOS[j];

							vector<double> jacobianos;

							jacobianos = transfCoordenadas(z,n);

							double detjac = jacobianos[0];

							double dz_i = dzpsi[gdli](z,n);
							double dn_i = dnpsi[gdli](z,n);

							double dfdx_i = dz_i * _J[0][0] + dn_i* _J[0][1];
							double dfdy_i = dz_i * _J[1][0] + dn_i* _J[1][1];

							double dz_j = dzpsi[gdlj](z,n);
							double dn_j = dnpsi[gdlj](z,n);

							double dfdx_j = dz_j * _J[0][0] + dn_j* _J[0][1];
							double dfdy_j = dz_j * _J[1][0] + dn_j* _J[1][1];

							KKUU += (C11 * dfdx_i * dfdx_j + C66 * dfdy_i * dfdy_j) * detjac * PESOS[j] * PESOS[i] * t;
							KKUV += (C12 * dfdx_i * dfdy_j + C66 * dfdy_i * dfdx_j) * detjac * PESOS[j] * PESOS[i] * t;
							KKVU += (C12 * dfdy_i * dfdx_j + C66 * dfdx_i * dfdy_j) * detjac * PESOS[j] * PESOS[i] * t;
							KKVV += (C11 * dfdy_i * dfdy_j + C66 * dfdx_i * dfdx_j) * detjac * PESOS[j] * PESOS[i] * t;
						}
					}
					KUU(gdli,gdlj)=KKUU;
					KUV(gdli,gdlj)=KKUV;
					KVU(gdli,gdlj)=KKVU;
					KVV(gdli,gdlj)=KKVV;
				}
			}
			MatrixXd K(2*N, 2*N);
			K<<KUU,KUV,KVU,KVV;
			guardarMatriz(rutaBase+"/Elemento"+to_string(nelemento)+"/KL_"+to_string(nelemento)+".csv",K);
		}
};