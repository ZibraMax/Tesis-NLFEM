#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include "legendre.h"


using namespace std;
int main(int argc, char const *argv[]) {
	cout << defaultfloat << setprecision(numeric_limits<double>::digits10);
	int n = 4;
	vector<double> puntos = darPuntos(n);
	vector<double> pesos = darPesos(n);
	for (int i = 0; i < n; ++i) {
		cout << puntos[i] << "," << pesos[i] << endl;
	}
	return 0;
}