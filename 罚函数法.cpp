#include <iostream>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

vector<double> advance(vector<double> x, vector<double> y, double h)
{
	vector<double> z(x.size());
	for (int i = 0; i < x.size(); i++)
	{
		z[i] = x[i] + y[i] * h;
	}
	return z;
}
double in_product(vector<double> x, vector<double> y)
{
	double r = 0;
	for (int i = 0; i < x.size(); i++)
	{
		r += x[i] * y[i];
	}
	return r;
}
double mod(vector<double> x)
{
	return pow(in_product(x, x), 0.5);
}
double max(double b)
{
	if (b<=0)
	{
		return 0;
	}
	else
	{
		return b;
	}
}
class result
{
public:
	vector<double> x;
	int n;
	bool k;
	double J;
	result(int nx, bool rk, double rJ);
	result(){};
	~result(){};
	void prt_result();
};

class penalty_function_method
{
public:
	vector<vector<double>> H, A;
	vector<double> f, b, lb, x;
	int nx, nc;
	double c, rho;
	penalty_function_method(double *h, double *i_f, double *a, double *i_b, double *x0,double ic,double rho0, int mx, int mc);
	~penalty_function_method(){};
	void prt_all();
	double val_J(vector<double> y);
	double val_phi(vector<double> y);
	vector<double> gradient();
	vector<double> optimal_next(vector<double> a, vector<double> dx, double h0);
	void next();
	void calculator();
};
penalty_function_method::penalty_function_method(double *h, double *i_f, double *a, double *i_b, double *x0, double ic, double rho0, int mx, int mc)
{
	nx = mx;
	nc = mc;
	rho = rho0;
	c = ic;
	H.resize(nx, vector<double>(nx));
	A.resize(nc, vector<double>(nx,0));
	f.resize(nx);
	b.resize(nc,0);
	lb.resize(nx, 0);
	x.resize(nx, 0);
	for (int i = 0; i < nx; i++)
	{
		f[i] = i_f[i];
		x[i] = x0[i];
		for (int j = 0; j < nx; j++)
		{
			H[i][j] = h[i*nx + j];
		}
	}
	for (int i = 0; i < nc; i++)
	{
		b[i] = i_b[i];
		for (int j = 0; j < nx; j++)
		{
			A[i][j] = a[i*nx + j];
		}
	}
}
void penalty_function_method::prt_all()
{
	cout << "H=" << endl;
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			cout << H[i][j] << '\t';
		}
		cout << endl;
	}
	cout << "f=" << endl;
	for (int i = 0; i < nx; i++)
	{
		cout << f[i] << '\t';
	}
	cout << endl << "A=" << endl;
	for (int i = 0; i < nc; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			cout << A[i][j] << '\t';
		}
		cout << endl;
	}
	cout << "b=" << endl;
	for (int i = 0; i < nc; i++)
	{
		cout << b[i] << '\t';
	}
	cout << endl << "x=" << endl;
	for (int i = 0; i < nx; i++)
	{
		cout << x[i] << '\t';
	}
	cout << endl;
}
double penalty_function_method::val_J(vector<double> y)
{
	double J = 0;
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			J = J + H[i][j] * y[i] * y[j]/2;
		}
		J += y[i] * f[i];
	}
	return J;
}
double penalty_function_method::val_phi(vector<double> y)
{
	double phi=0, temp;
	for (int i = 0; i < nc; i++)
	{
		if (i < nx)
		{
			for (int j = 0; j < nx; j++)
			{
				phi = phi + H[i][j] * y[i] * y[j];
			}
			phi += y[i] * f[i];
		}
		temp = 0;
		for (int j = 0; j < nx; j++)
		{
			temp += A[i][j] * y[j];
		}
		phi += rho*pow(max(temp - b[i]), 2);
	}
	return phi;
}
vector<double> penalty_function_method::gradient()
{
	vector<double> DJ(f);
	double temp;
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			if (j<nx)
			{
				DJ[i] = DJ[i] + H[i][j] * x[j];
			}
			temp = 0;
			for (int k = 0; k < nx; k++)
			{
				temp += A[j][k] * x[k];
			}
			DJ[i] += 2 * max(temp - b[j])*A[i][i] * x[i];
		}
		DJ[i] = -DJ[i];
	}
	return DJ;
}
vector<double> penalty_function_method::optimal_next(vector<double> a, vector<double> dx, double h0)
{
	double h = h0;
	a = advance(a, dx, h);
	vector<double> b = advance(a, dx, h);
	for (int j = 0;; j++)
	{
		for (int i = 0; val_phi(a)>val_phi(b); i++)
		{
			h = h * 2;
			a = b;
			b = advance(a, dx, h);
		}
		if (h<0.01)
		{
			break;
		}
		h = h0;
		b = advance(a, dx, h);
	}
	return a;
}
void penalty_function_method::next()
{
	vector<double> dx = gradient(), next_x = optimal_next(x, dx, 0.001);
	x = next_x;
	rho = rho*c;
}
void penalty_function_method::calculator()
{
	double l = 1;
	for (int i = 0;l>pow(0.1, 4); i++)
	{
		next();
		l = mod(gradient());
	}

	cout << val_J(x) << endl;
}

int main()
{
	double h[4] = { 6, -4, -4, 4 }, f[2] = { 3, -4 }, a[8] = { 2, 1, -1, 2, -1, 0, 0, -1 }, b[4] = { 4, 4, 0, 0 }, x[4] = { 1, 1 };
	penalty_function_method C(h, f, a, b, x,10,100, 2, 4);
	C.calculator();
	system("pause");
	return 0;
}