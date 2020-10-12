#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Serendipity.h"
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

double t = 0.5;
double E = 2.1*pow(10.0,6);
double v = 0.2;

double l = 0.1;
double C11 = E/(1.0-v*v);
double C12 = v*E/(1.0-v*v);
double C66 = E/2/(1.0+v);
const double PI_calc  = 3.141592653589793238463;
double L0 = (1.0)/(2*PI_calc*l*l*t);



double PUNTOS[] = {-sqrt(3.0/5.0),0,sqrt(3.0/5.0)};
int NG = sizeof(PUNTOS)/sizeof(PUNTOS[0]);
double PESOS[] = {5.0/9.0,8.0/9.0,5.0/9.0};

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

vector<Serendipity> leerTexto(string filename) {
	vector<Serendipity> ELEMENTOS;
	string line;
	ifstream myfile (filename);
	if (myfile.is_open()) {
		getline(myfile,line);
		vector<string> linea = split(line," ");
		int NUMERO_GDL = stoi(linea[0]);
		int NUMERO_ELEMENTOS = stoi(linea[1]);
		int NUMERO_NODOS_ELEMENTO = stoi(linea[2]);
		double GDL[NUMERO_GDL][2];
		vector<vector<int>> elementos;
		vector<int> ENL[NUMERO_ELEMENTOS];
		for (int i = 0; i < NUMERO_GDL; ++i) {
			getline(myfile,line);
			vector<string> linea = split(line," ");
			GDL[i][0] = stod(linea[0]);
			GDL[i][1] = stod(linea[1]);
		}
		for (int i = 0; i < NUMERO_ELEMENTOS; ++i) {
			getline(myfile,line);
			vector<string> linea = split(line," ");
			vector<int> elemento;
			for (int j = 0; j < NUMERO_NODOS_ELEMENTO; ++j) {
				elemento.push_back(stoi(linea[j])-1);
			}
			elementos.push_back(elemento);
		}
		for (int i = 0; i < NUMERO_ELEMENTOS; ++i) {
			getline(myfile,line);
			vector<string> linea = split(line,"\t");
			int L = stoi(linea[0]);
			vector<int> nolocales;
			for (int j = 1; j <= L; ++j) {
				nolocales.push_back(stoi(linea[j])-1);
			}
			ELEMENTOS.push_back(Serendipity(elementos[i],nolocales,GDL));
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	return ELEMENTOS;
}
void writeToCSVfile(string name, MatrixXd matrix) {
	cout.precision(15);
	const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
	ofstream Archivo(name.c_str());
    Archivo << matrix.format(CSVFormat);
	Archivo.close();
}
double atenuacion(double x,double y,double xnl,double ynl) {
	double distancia = sqrt(pow((x-xnl),2)+pow((y-ynl),2));
	return L0*exp(-distancia/l);
}
int main () {
	cout<<PI_calc<<endl;
	vector<Serendipity> ELEMENTOS;
	ELEMENTOS = leerTexto("input.txt");
	int conti = 0;
	string Ruta = "./MATRICES"; 
    int n = Ruta.length(); 
    char RutaChar[n + 1]; 
    strcpy(RutaChar, Ruta.c_str());
	mkdir(RutaChar);

	for (Serendipity e: ELEMENTOS) {
		int contj = 0;
		conti++;
		string Ruta = "./MATRICES/Elemento"+to_string(conti); 
	    int n = Ruta.length(); 
	    char RutaChar[n + 1]; 
	    strcpy(RutaChar, Ruta.c_str());
		mkdir(RutaChar);
		for (int indiceNoLocal: e.nolocales) {
			contj++;
			Serendipity enl = ELEMENTOS[indiceNoLocal];
			int N = e.gdl.size();
			int NL = enl.gdl.size();
			MatrixXd KUU(N,NL);
			MatrixXd KUV(N,NL);
			MatrixXd KVU(N,NL);
			MatrixXd KVV(N,NL);
			for (int gdli = 0; gdli < e.gdl.size(); ++gdli) {
				for (int gdlj = 0; gdlj < enl.gdl.size(); ++gdlj) {
					double KKUU = 0.0;
					double KKUV = 0.0;
					double KKVU = 0.0;
					double KKVV = 0.0;
					for (int i = 0; i < NG; ++i) {
						double z = PUNTOS[i];
						for (int j = 0; j < NG; ++j) {
							double n = PUNTOS[j];
							double x = e.TX(z,n);
							double y = e.TY(z,n);
							vector<double> jacobianos;
							jacobianos = e.transfCoordenadas(z,n);
							double detjac = jacobianos[0];
							vector<vector<double>> jacobiano_ = e._J;

                            vector<double> dzpsis = e.dzpsis(z,n);
                            vector<double> dnpsis = e.dnpsis(z,n);

                            double dz_i = dzpsis[gdli];
                            double dn_i = dnpsis[gdli];

                            double dfdx_i = dz_i * jacobiano_[0][0] + dn_i* jacobiano_[0][1];
                            double dfdy_i = dz_i * jacobiano_[1][0] + dn_i* jacobiano_[1][1];

							for (int i_nl = 0; i_nl < NG; ++i_nl) {
								double znl = PUNTOS[i_nl];
								for (int j_nl = 0; j_nl < NG; ++j_nl) {
									double nnl = PUNTOS[j_nl];
									double xnl = enl.TX(znl,nnl);
									double ynl = enl.TY(znl,nnl);
									vector<double> jacobianosnl;

									jacobianosnl = enl.transfCoordenadas(znl,nnl);
									vector<vector<double>> jacobiano_nl = enl._J;
									double detjacnl = jacobianosnl[0];

		                            vector<double> dzpsisnl = enl.dzpsis(znl,nnl);
		                            vector<double> dnpsisnl = enl.dnpsis(znl,nnl);

		                            double dznl_j = dzpsisnl[gdlj];
		                            double dnnl_j = dnpsisnl[gdlj];

		                            double dfdxnl_j = dznl_j * jacobiano_nl[0][0] + dnnl_j* jacobiano_nl[0][1];
		                            double dfdynl_j = dznl_j * jacobiano_nl[1][0] + dnnl_j* jacobiano_nl[1][1];
									double AZN = atenuacion(x,y,xnl,ynl);
									KKUU += AZN * (C11 * dfdx_i * dfdxnl_j + C66 * dfdy_i * dfdynl_j) * detjac * detjacnl * PESOS[j_nl] * PESOS[i_nl] * PESOS[j] * PESOS[i];
									KKUV += AZN * (C12 * dfdx_i * dfdynl_j + C66 * dfdy_i * dfdxnl_j) * detjac * detjacnl * PESOS[j_nl] * PESOS[i_nl] * PESOS[j] * PESOS[i];
									KKVU += AZN * (C12 * dfdy_i * dfdxnl_j + C66 * dfdx_i * dfdynl_j) * detjac * detjacnl * PESOS[j_nl] * PESOS[i_nl] * PESOS[j] * PESOS[i];
									KKVV += AZN * (C11 * dfdy_i * dfdynl_j + C66 * dfdx_i * dfdxnl_j) * detjac * detjacnl * PESOS[j_nl] * PESOS[i_nl] * PESOS[j] * PESOS[i];
								}
							}
						}
					}
					KUU(gdli,gdlj)=KKUU;
					KUV(gdli,gdlj)=KKUV;
					KVU(gdli,gdlj)=KKVU;
					KVV(gdli,gdlj)=KKVV;
				}
			}
			MatrixXd K(2*N, 2*NL);
			K<<KUU,KUV,KVU,KVV;
			writeToCSVfile("./MATRICES/Elemento"+to_string(conti)+"/KNL"+to_string(conti)+"_"+to_string(contj)+".csv",K);
		}
		cout<<"Local: "<<conti<<endl;
	}
	return 0;
}