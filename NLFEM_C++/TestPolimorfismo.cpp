// NLFEM - Integrador
// Arturo Rodriguez - da.rodriguezh@uniandes.edu.co
// Basado en "Nonlocal integral elasticity: 2D finite element based solutions" (Pisano, 2009) y en el enmallado de Fernando Ramírez
// Programa usado SOLAMENTE para integrar.

// Librerias Requeridas:
// Eigen (paquete de algebra lineal para C++)


// Inputs:

/* Argumentos de entrada (por consola separados por espacio): 
1. Archivo de texto que contenga las siguientes filas:
	1 fila: Numero de grados de libertad, Numero de elementos
	... : Coordenadas (2D) de los grados de libertad
	... : Grados de libertad de cada uno de los elementos
	... : Elementos no locales por elemento
2. Numero de puntos de gauss para integrar (int)
3. Modulo de Young (double)
4. Coeficiente de Poisson (double)
5. grosor (double)
6. longitud interna (parámetro no-local) (double)
7. Nombre de la carpeta donde se quiere guardar (Se creara de no existir) (string)
8. Funcion de atenuación a usar. Entero entre 1 y 4 asi:
	1: Función Biexponencial
	2: Función Cónica (lineal)
	3: Función campana (Cuadrática)
	4: Función custom (biexponencial modificada)
*/

// Outputs:

/* El programa generará (y sobreescribirá de ser necesario) los resultados en la siguiente estructura

	index.cpp
	Serendipity.h
	MATRICES
	--Elementoi
	----KL_i.csv
	----KNLS.csv

	El archivo KL_i.csv contiene la matriz de rigidez local del elemento i, esta se guarda en forma de matriz (16x16)
	El archivo KNLS.csv contiene M filas donde cada j-ésima fila corresponde a la matriz de rigidez no local
	del elemento i con el elemento j. Esta matriz se guarda por columnas. Es decir, la posicion 18 corresponde a la fila 2 columna 2 de la matriz.
	Para poder cargar cada una de las matrices no locales (en python) se puede usar el método reshape de numpy y
	posteriormente hallar la transpuesta, en numpy asi: (elemento i) -> KNLS[j] = np.loadtxt('MATRICES/KNLS.csv')[j].reshape([16,16]).T
*/

/*
	Los argumentos de entrada permiten llamar al programa desde python, lo cual puede ser util.
*/

/*
	En ningun punto de este programa se usan variables especificamente relacionadas con el numero de nodos por elemento
	Esto hace que usando otras clases (como Serendipity.h) se pueden implementar elementos con 4 nodos por elemento,
	asi como elementos triangulares de 3 y 6 nodos.

	En algún punto se tendrá que definir la interface Elemento con el polimorfismo para cada una de estas clases


	Adicionalmente se tiene un contador de tiempo asociado al tiempo restante de integración.
	Este tiempo restante es solamente confiable cuando se complenta un 20% de los elementos,
	esto se debe a que los primeros elementos del enmallado generalmente tienen una menor cantidad de elementos no locales
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Serendipity.h"
#include "TriangularC.h"
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include "Eigen/Dense"
#include <chrono> 
#include <iomanip>
#include <limits>
#include "legendre.h"


using namespace std::chrono;
using namespace Eigen;
using namespace std;

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

vector<Elemento*> leerTexto(string filename,double t,vector<double> PUNTOS, vector<double> PESOS) {
	vector<Elemento*> ELEMENTOS;
	string line;
	ifstream myfile (filename);
	if (myfile.is_open()) {
		getline(myfile,line);
		vector<string> linea = split(line,"\t");
		int NUMERO_GDL = stoi(linea[0]);
		int NUMERO_ELEMENTOS = stoi(linea[1]);
		double GDL[NUMERO_GDL][2];
		vector<vector<int>> elementos;
		vector<int> ENL[NUMERO_ELEMENTOS];
		for (int i = 0; i < NUMERO_GDL; ++i) {
			getline(myfile,line);
			vector<string> linea = split(line,"\t");
			GDL[i][0] = stod(linea[0]);
			GDL[i][1] = stod(linea[1]);
		}
		for (int i = 0; i < NUMERO_ELEMENTOS; ++i) {
			getline(myfile,line);
			vector<string> linea = split(line,"\t");
			vector<int> elemento;
			for (int j = 0; j < linea.size(); ++j) {
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
			if (elementos[i].size()==8) {
				ELEMENTOS.push_back(new Serendipity(elementos[i],nolocales,GDL,PUNTOS,PESOS,t));
			} else if (elementos[i].size()==6) {
				ELEMENTOS.push_back(new TriangularC(elementos[i],nolocales,GDL,PUNTOS,PESOS,t));
			} else {
				cout<<"Error al cargar el elemento " <<i<<endl;
			}
		}
		myfile.close();
		cout<<"Se termino de cargar la geometria"<<endl;
	}
	else cout << "No se pudo abrir el archivo o no se encuentra";
	return ELEMENTOS;
}

//Función de atenuación que se usará
double atenuacion(double x,double y,double xnl,double ynl,double L0,double l) {
	double distancia = sqrt(pow((x-xnl),2)+pow((y-ynl),2));
	return L0*exp(-distancia/l);
}
double atenuacion_biexponencial_modificada(double x,double y,double xnl,double ynl,double L0,double l) {
	double distancia = sqrt(pow((x-xnl),2)+pow((y-ynl),2));
	return L0*distancia/l*exp(-distancia/l);
}
double atenuacion_lineal(double x,double y,double xnl,double ynl,double L0,double LR) {
	double r = sqrt(pow((x-xnl),2)+pow((y-ynl),2));
	return (r/LR) > 1 ? 0 : L0*(1-(r/LR));
}
double atenuacion_cuadratica(double x,double y,double xnl,double ynl,double L0,double LR) {
	double r = sqrt(pow((x-xnl),2)+pow((y-ynl),2));
	return (pow(r,2)/pow(LR,2)) > 1 ? 0 : L0*(1-(pow(r,2)/pow(LR,2)));
}
int main (int argc, char const *argv[]) {

	cout<<"NLFEM-C++ Arturo Rodriguez & Fernando Ramirez"<<endl;
	cout << "================================" <<endl;
	int NG = stoi(argv[2]);
	cout << "Encontrando puntos y pesos de Gauss" <<endl;
	vector<double> PUNTOS = darPuntos(NG);
	vector<double> PESOS = darPesos(NG);
	cout << "Puntos y pesos encontrados (" << NG << "):" <<endl;
	for (int i = 0; i < NG; ++i) {
		cout << "Z"<<i<<"=" << PUNTOS[i] << ", W"<<i<<"=" << PESOS[i] << endl;
	}
	cout << "--------------------------------" <<endl;
	cout << "Calculando con los siguientes parametros:" <<endl;
	double E = stold(argv[3]);
	double v = stold(argv[4]);
	double t = stold(argv[5]);
	double l = stold(argv[6]);
	cout<< "E=" << E << ","<< "v="<< v << ","<< "t=" << t << "," << "l=" << l << endl;
	cout << "--------------------------------" <<endl;
	cout << "Archivo de enmallado: " << argv[1] <<endl;
	double C11 = E/(1.0-v*v);
	double C12 = v*E/(1.0-v*v);
	double C66 = E/2/(1.0+v);
	const double PI_calc  = 3.141592653589793238463;

	vector<Elemento*> ELEMENTOS;
	cout << "Leyendo archivo" << endl;
	ELEMENTOS = leerTexto(argv[1],t,PUNTOS,PESOS);
	cout << "--------------------------------" <<endl;
	cout << "Comenzando la integracion" << endl;
	int NUMERO_ELEMENTOS = ELEMENTOS.size();
	int conti = 0;
	string rutaParcial = argv[7];
	string Ruta = "./"+rutaParcial; 
	int n = Ruta.length();
	char RutaChar[n + 1];
	strcpy(RutaChar, Ruta.c_str());
	mkdir(RutaChar);
	double promedio = 0;
	double tiempo = 0;
	for (Elemento *e: ELEMENTOS) {
		auto start = high_resolution_clock::now();
		int contj = 0;
		conti++;
		string Ruta = "./"+rutaParcial+"/Elemento"+to_string(conti); 
		int n = Ruta.length();
		char RutaChar[n + 1];
		strcpy(RutaChar, Ruta.c_str());
		mkdir(RutaChar);
		e->matrizLocal(C11,C12,C66,conti,"./"+rutaParcial);
		vector<vector<double>> matricesNoLocales;
		for (int indiceNoLocal: e->nolocales) {
			Elemento *enl = ELEMENTOS[indiceNoLocal];
			contj++;
			int N = e->gdl.size();
			int NL = enl->gdl.size();
			MatrixXd KUU(N,NL);
			MatrixXd KUV(N,NL);
			MatrixXd KVU(N,NL);
			MatrixXd KVV(N,NL);
			for (int gdli = 0; gdli < e->gdl.size(); ++gdli) {
				for (int gdlj = 0; gdlj < enl->gdl.size(); ++gdlj) {
					double KKUU = 0.0;
					double KKUV = 0.0;
					double KKVU = 0.0;
					double KKVV = 0.0;
					for (int i = 0; i < e->PESOS.size(); ++i) {
						double z = e->PUNTOS[i][0];
						double n = e->PUNTOS[i][1];
						double x = e->TX(z,n);
						double y = e->TY(z,n);
						vector<double> jacobianos;
						jacobianos = e->transfCoordenadas(z,n);
						double detjac = jacobianos[0];
						vector<vector<double>> jacobiano_ = e->_J;
						double dz_i = e->dzpsi[gdli](z,n);
						double dn_i = e->dnpsi[gdli](z,n);
						double dfdx_i = dz_i * jacobiano_[0][0] + dn_i* jacobiano_[0][1];
						double dfdy_i = dz_i * jacobiano_[1][0] + dn_i* jacobiano_[1][1];
						for (int i_nl = 0; i_nl < enl->PESOS.size(); ++i_nl) {
							double znl = enl->PUNTOS[i_nl][0];
							double nnl = enl->PUNTOS[i_nl][1];
							double xnl = enl->TX(znl,nnl);
							double ynl = enl->TY(znl,nnl);
							vector<double> jacobianosnl;

							jacobianosnl = enl->transfCoordenadas(znl,nnl);
							vector<vector<double>> jacobiano_nl = enl->_J;
							double detjacnl = jacobianosnl[0];

							double dznl_j = enl->dzpsi[gdlj](znl,nnl);
							double dnnl_j = enl->dnpsi[gdlj](znl,nnl);

							double dfdxnl_j = dznl_j * jacobiano_nl[0][0] + dnnl_j* jacobiano_nl[0][1];
							double dfdynl_j = dznl_j * jacobiano_nl[1][0] + dnnl_j* jacobiano_nl[1][1];
							double AZN = 0;
							if (stoi(argv[8])==1) {
								double L0 = 1.0/2/PI_calc/l/l/t;
								AZN = atenuacion(x,y,xnl,ynl,L0,l);
							} else if (stoi(argv[8])==2) {
								double LR =6*l;
								double L0 = 3.0/2.0/PI_calc/LR/LR/t;
								AZN = atenuacion_lineal(x,y,xnl,ynl,L0,LR);
							} else if (stoi(argv[8])==3) {
								double LR =6*l;
								double L0 = 1.0/PI_calc/LR/LR/t;
								AZN = atenuacion_cuadratica(x,y,xnl,ynl,L0,LR);
							} else if (stoi(argv[8])==4) {
								double L0 = 1.0/4.0/PI_calc/l/l/t;
								AZN = atenuacion_biexponencial_modificada(x,y,xnl,ynl,L0,l);
							} else {
								cout<<"ERROR EN ELECCION DE FUNCION DE ATENUACION";
							}
							KKUU += t * t * AZN * (C11 * dfdx_i * dfdxnl_j + C66 * dfdy_i * dfdynl_j) * detjac * detjacnl * enl->PESOS[i_nl] * e->PESOS[i];
							KKUV += t * t * AZN * (C12 * dfdx_i * dfdynl_j + C66 * dfdy_i * dfdxnl_j) * detjac * detjacnl * enl->PESOS[i_nl] * e->PESOS[i];
							KKVU += t * t * AZN * (C12 * dfdy_i * dfdxnl_j + C66 * dfdx_i * dfdynl_j) * detjac * detjacnl * enl->PESOS[i_nl] * e->PESOS[i];
							KKVV += t * t * AZN * (C11 * dfdy_i * dfdynl_j + C66 * dfdx_i * dfdxnl_j) * detjac * detjacnl * enl->PESOS[i_nl] * e->PESOS[i];
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
			Map<const VectorXd> v1(K.data(), K.size());
			vector<double> v2;
			v2.resize(v1.size());
			VectorXd::Map(&v2[0], v1.size()) = v1;
			matricesNoLocales.push_back(v2);
		}
		ofstream myfile("./"+rutaParcial+"/Elemento"+to_string(conti)+"/KNLS.csv");
		myfile << defaultfloat << setprecision(numeric_limits<double>::digits10);
		for (int i = 0; i < matricesNoLocales.size(); ++i) {
			for (int j = 0; j < matricesNoLocales[i].size(); ++j) {
				if (j == matricesNoLocales[i].size()-1) {
					myfile << matricesNoLocales[i][j];
				} else {
					myfile << matricesNoLocales[i][j] << ",";
				}
			}
			myfile << endl;
		}
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		tiempo+=duration.count();
		promedio =  (promedio*(conti-1)+duration.count())/conti;
		double porcentaje = 100.0*conti/NUMERO_ELEMENTOS;
		cout<<"Local: "<<conti<<" - Tiempo: "<< duration.count()<<"ms - "<< porcentaje << "%"<< " - ETA: " << ceil(promedio*NUMERO_ELEMENTOS-promedio*conti)/1000/60 << " minutos"<<endl;
	}
	cout<<"Se ha terminado de calcular. Tiempo trascurrido: "<<tiempo/1000/60<<" min"<<endl;
	return 0;
}