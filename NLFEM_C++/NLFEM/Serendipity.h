#include <vector>
#include <iostream>
using namespace std;

class Serendipity {
	public:
		vector<vector<double>> coords;
		vector<int> gdl;
		vector<int> nolocales;
		vector<double> X;
		vector<double> Y;
		vector<vector<double>> J;
		vector<vector<double>> _J;
		Serendipity(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2]) {
			gdl = pgdls;
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
				(double) 1/4*(1-z)*(1-n)*(-1-z-n),
				(double) 1/4*(1+z)*(1-n)*(-1+z-n),
				(double) 1/4*(1+z)*(1+n)*(-1+z+n),
				(double) 1/4*(1-z)*(1+n)*(-1-z+n),
				(double) 1/2*(1-z*z)*(1-n),
				(double) 1/2*(1+z)*(1-n*n),
				(double) 1/2*(1-z*z)*(1+n),
				(double) 1/2*(1-z)*(1-n*n)
				};
			return foo;
		}
		vector<double> dzpsis(double z, double n) {
			vector<double> foo{
				(double) -1/4*(n-1)*(2*z+n),
				(double) -1/4*(n-1)*(2*z-n),
				(double) 1/4*(n+1)*(2*z+n),
				(double) 1/4*(n+1)*(2*z-n),
				(double) (n-1)*z,
				(double) -1/2*(n*n-1),
				(double) -(n+1)*z,
				(double) 1/2*(n*n-1)
				};
			return foo;
		}
		vector<double> dnpsis(double z, double n) {
			vector<double> foo{
				(double) -1/4*(z-1)*(2*n+z),
				(double) 1/4*(z+1)*(2*n-z),
				(double) 1/4*(z+1)*(2*n+z),
				(double) -1/4*(z-1)*(2*n-z),
				(double) 1/2*(z*z-1),
				(double) -n*(z+1),
				(double) -1/2*(z*z-1),
				(double) n*(z-1)
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
};