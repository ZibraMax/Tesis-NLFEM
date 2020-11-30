#include <vector>
#include <cmath>

using namespace std;

const double PI_calc  = 3.141592653589793238463;

double p(double x,int n) {
	if (n==0) {
		return 1;
	} else if(n==1) {
		return x;
	} else {
		return 1.0/n*((2.0*n-1.0)*x*p(x,n-1)-(n-1.0)*p(x,n-2));
	}
}
double dpdx(double x,int n) {
	return n/(x*x-1.0)*(x*p(x,n)-p(x,n-1));
}

vector<double> darPuntos(int n) {
	vector<double> puntos;
	double tol = 1.0*pow(10.0,-16.0);
	for (int i = 1; i <= n; ++i) {
		double x = cos(PI_calc * (i - 0.25) / (n + 0.5));
		double error = 1;
		while (abs(error) > tol) {
			error = p(x,n)/dpdx(x,n);
			x -= error;
		}
		puntos.push_back(x);
	}
	return puntos;
}

vector<double> darPesos(int n) {
	vector<double> puntos = darPuntos(n);
	vector<double> pesos;
	for (int i = 0; i < puntos.size(); ++i) {
		double deriv = dpdx(puntos[i],n);
		double w = 2.0/(1.0-puntos[i]*puntos[i])/deriv/deriv;
		pesos.push_back(w);		
	}
	return pesos;
}