#include <vector>
#include <iostream>
using namespace std;

class Serendipity {
	public:
		vector<vector<double>> vertices;
		vector<int> gdl;
		vector<int> nolocales;
		Serendipity(vector<int> pgdls, vector<int> pnolocales,double COORDENADAS[][2]) {
			gdl = pgdls;
			nolocales = pnolocales;
			vertices = darCoordenadas(COORDENADAS,gdl);
		}
		vector<vector<double>> darCoordenadas(double COORDENADAS[][2], vector<int> gdls) {
			vector<vector<double>> v;
			for (int i = 0; i < gdls.size(); ++i) {
				vector<double> v1{COORDENADAS[gdls[i]][0],COORDENADAS[gdls[i]][1]};
				v.push_back(v1);
			}
			return v;
		}
		void imprimirse() {
			for (vector<double> vertice: vertices) {
				for (double coordenada: vertice) {
					cout<< coordenada << "\t";
				}
				cout << endl;
			}
		}
};