#include <iostream>
#include <vector>
#include <string>
#include <math.h>
using namespace std;

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
result::result(int nx, bool rk,double rJ)
{
	n = nx;
	k = rk;
	x.resize(nx);
	J = rJ;
}
void result::prt_result()
{
	if (k)
	{
		cout << "最优解为：" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << "x"<<i+1<<"="<<x[i] << '\t';
		}
		cout << endl;
		cout << "Jmin=" << J << endl;
	}
	else
	{
		cout << "不存在最优解" << endl;
	}
}
class simplex_method
{
public:
	vector<vector<double>> S;
	vector<int> unbase;
	int nx, nc;
	simplex_method(double *f, double *a, double *b, int mx, int mc);
	simplex_method(){};
	~simplex_method(){};
	void prt_S();
	result put_result(double J=0);
	vector<double> cal_J();
	int enter_base(vector<double> J);
	int leave_base(int e);
	void guass_eli(int i,int j);
	void swap_row(int i, int j);
	void ex_base(int e, int l);
	void calcultor(result &r);
};
simplex_method::simplex_method(double *f, double *a, double *b, int mx, int mc)
{
	nx = mx;
	nc = mc;
	S.resize(nc + 1, vector<double>(nx + nc + 2, 0));
	unbase.resize(nx);
	for (int i = 0; i < nx; i++)
	{
		unbase[i] = i + 1;
	}
	for (int i = 0; i < nx; i++)
	{
		S[0][i + 1] = f[i];
	}
	for (int i = 1; i < nc+1; i++)
	{
		S[i][0] = nx + i;
		S[i][nx + nc + 1] = b[i-1];
		S[i][i+nx] = 1;
		for (int j = 1; j < nx+1; j++)
		{
			S[i][j] = a[(i-1)*nx + j-1];
		}
	}
}
void simplex_method::prt_S()
{
	for (int i = 0; i < nc+1; i++)
	{
		for (int j = 0; j < nx+nc+2; j++)
		{
			cout << S[i][j] << '\t';
		}
		cout << endl;
	}
}
result simplex_method::put_result(double J)
{
	result r(nx,S[0][0],J);
	for (int i = 0; i < nx; i++)
	{
		r.x[i] = S[i + 1][nx + nc + 1];
	}
	return r;
}
vector<double> simplex_method::cal_J()
{
	vector<double> J(nx);
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			J[i] = J[i]+S[j + 1][unbase[i]] * S[0][S[j + 1][0]];
		}
	}
	for (int i = 0; i < nx; i++)
	{
		cout << unbase[i] << '\t';
	}
	cout << endl;
	for (int i = 0; i < nx; i++)
	{
		cout << J[i] << '\t';
	}
	cout << endl << endl;
	return J;
}
int simplex_method::enter_base(vector<double> J)
{
	double sigma;
	int k = 1;
	sigma = S[0][unbase[0]] - J[0];
	for (int i = 1; i < nx; i++)
	{
		if (sigma > S[0][unbase[i]] - J[i])
		{
			sigma = S[0][unbase[i]] - J[i];
			k = unbase[i];
		}
	}
	if (sigma>=0)
	{
		return 0;
	}
	return k;
}
int simplex_method::leave_base(int e)
{
	bool b = 0;
	for (int i = 1; i < nc+1; i++)
	{
		if (S[i][e]>0) b = 1;
	}
	if (!b)
	{
		return -1;
	}
	double temp;
	int l=1;
	temp = S[1][nx + nc + 1] / S[1][e];
	for (int i = 2; i < nc+1; i++)
	{
		if (temp > S[i][nx + nc + 1] / S[i][e])
		{
			temp = S[i][nx + nc + 1] / S[i][e];
			l = i;
		}
	}
	if (e!=0)
	{
		for (int i = 0; i < nx; i++)
		{
			if (unbase[i] == e)
			{
				unbase[i] = S[l][0];
			}
		}
		S[l][0] = e;
	}
	return l;
}
void simplex_method::guass_eli(int i,int j)
{
	double temp_ij = S[i][j];
	for (int k = 1; k < nx+nc+2; k++)
	{
		S[i][k] = S[i][k] / temp_ij;
	}
	for (int k = 1; k < nc+1; k++)
	{
		double temp = S[k][j];
		for (int l = 1; l < nx+nc+2; l++)
		{
			if (k != i)
			{
				S[k][l] = S[k][l] - S[i][l] * temp;
			}
		}
	}
}
void simplex_method::ex_base(int e, int l)
{
	guass_eli(l, e);
}
void simplex_method::swap_row(int i, int j)
{
	double temp;
	for (int k = 0; k < nx+nc+2; k++)
	{
		temp = S[i][k];
		S[i][k] = S[j][k];
		S[j][k] = temp;
	}
}
void simplex_method::calcultor(result &r)
{
	vector<double> J(nx);
	double f=0;
	int e=1, l=1;
	for (int i = 0; ; i++)
	{
		J = cal_J();
		e = enter_base(J);
		l = leave_base(e);
		if (e == 0 || l == -1) break;
		ex_base(e, l);
		prt_S();
		cout << endl;
	}
	if (e==0)
	{
		S[0][0] = 1;
	}
	for (int i = 1; i < nx+1; i++)
	{
		for (int j = 1; j < nc+1; j++)
		{
			if (S[j][0] == i) swap_row(i, j);
		}
	}
	for (int i = 0; i < nx; i++)
	{
		f = f + S[0][i + 1] * S[i + 1][nx + nc + 1];
	}
	prt_S();
	r = put_result(f);
}
int main()
{
	double f[2] = { -3, -4 }, a[6] = { 2, 1, 1, 3 }, b[3] = { 40, 30 };
	simplex_method H(f, a, b, 2, 2);
	H.prt_S();
	cout << endl;
	result r;
	H.calcultor(r);
	r.prt_result();
	system("pause");
	return 0;
}