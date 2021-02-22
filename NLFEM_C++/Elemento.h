#ifndef ELEMENTO_H
#define ELEMENTO_H
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

class Elemento {
	public:
		vector<vector<double>> coords;
		vector<int> gdl;
		vector<int> nolocales;
		vector<double> X;
		vector<double> Y;
		vector<vector<double>> J;
		vector<vector<double>> _J;
		vector<vector<double>> PUNTOS;
		vector<double> PESOS;
		double t;

		vector<function<double(double,double)>> psi;
		vector<function<double(double,double)>> dzpsi;
		vector<function<double(double,double)>> dnpsi;

		Elemento(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2],vector<double> ppgauss,vector<double> pwgauss,double pt) {
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
			for (int i = 0; i < psi.size(); ++i) {
				res+=X[i]*psi[i](z,n);
			}
			return res;
		}
		double TY(double z,double n) { 
			double res = 0;
			for (int i = 0; i < psi.size(); ++i) {
				res+=Y[i]*psi[i](z,n);
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
		vector<double> transfCoordenadas(double z,double n) {
			double dxdz = 0.0;
			double dydz = 0.0;
			double dxdn = 0.0;
			double dydn = 0.0;
			for (int i = 0; i < coords.size(); ++i) {
				dxdz += X[i]*dzpsi[i](z,n);
				dydz += Y[i]*dzpsi[i](z,n);
				dxdn += X[i]*dnpsi[i](z,n);
				dydn += Y[i]*dnpsi[i](z,n);
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
						double z = PUNTOS[i][0];
						double n = PUNTOS[i][1];

						vector<double> jacobianos;

						jacobianos = transfCoordenadas(z,n);

						// cout<<jacobianos[0]<<","<<jacobianos[1]<<endl;

						double detjac = jacobianos[0];

						double dz_i = dzpsi[gdli](z,n);
						double dn_i = dnpsi[gdli](z,n);

						double dfdx_i = dz_i * _J[0][0] + dn_i* _J[0][1];
						double dfdy_i = dz_i * _J[1][0] + dn_i* _J[1][1];

						double dz_j = dzpsi[gdlj](z,n);
						double dn_j = dnpsi[gdlj](z,n);

						double dfdx_j = dz_j * _J[0][0] + dn_j* _J[0][1];
						double dfdy_j = dz_j * _J[1][0] + dn_j* _J[1][1];

						KKUU += (C11 * dfdx_i * dfdx_j + C66 * dfdy_i * dfdy_j) * detjac * PESOS[i] * t;
						KKUV += (C12 * dfdx_i * dfdy_j + C66 * dfdy_i * dfdx_j) * detjac * PESOS[i] * t;
						KKVU += (C12 * dfdy_i * dfdx_j + C66 * dfdx_i * dfdy_j) * detjac * PESOS[i] * t;
						KKVV += (C11 * dfdy_i * dfdy_j + C66 * dfdx_i * dfdx_j) * detjac * PESOS[i] * t;
					}

					KUU(gdli,gdlj)=KKUU;
					KUV(gdli,gdlj)=KKUV;
					KVU(gdli,gdlj)=KKVU;
					KVV(gdli,gdlj)=KKVV;
				}
			}
			MatrixXd K(2*N, 2*N);
			K<<KUU,KUV,KVU,KVV;
			cout.precision(15);
			const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
			ofstream Archivo((rutaBase+"/Elemento"+to_string(nelemento)+"/KL_"+to_string(nelemento)+".csv").c_str());
		    Archivo << K.format(CSVFormat);
			Archivo.close();
		}
};
#endif