#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Serendipity.h"
#include <cstdlib>
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
				elemento.push_back(stoi(linea[j]));
			}
			elementos.push_back(elemento);
		}
		for (int i = 0; i < NUMERO_ELEMENTOS; ++i) {
			getline(myfile,line);
			vector<string> linea = split(line,"\t");
			int L = stoi(linea[0]);
			vector<int> nolocales;
			for (int j = 0; j <= L; ++j) {
				nolocales.push_back(stoi(linea[j]));
			}
			ELEMENTOS.push_back(Serendipity(elementos[i],nolocales,GDL));
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	return ELEMENTOS;
}

int main () {
	vector<Serendipity> ELEMENTOS;
	ELEMENTOS = leerTexto("input.txt");
	ofstream Archivo("elementos.csv");
	for (Serendipity elemento: ELEMENTOS) {
		for (vector<double> vertice: elemento.vertices) {
			for (double coordenada: vertice) {
				Archivo << coordenada << ",";
			}
		}
		Archivo << endl;
	}
	Archivo.close();
	return 0;
}