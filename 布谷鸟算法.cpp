#include<iostream>
#include<math.h>
#include<vector>
#include "stdlib.h"
#include <ctime>
#include <fstream>
#include <string>
using namespace std;
#define pi 3.141592653
#define lamd 1.5
#define alpha 1
#define af 0.5
double f(vector<double> x)
{
	double y=0,l=1;
	for (unsigned i = 0; i < x.size(); i++)
	{
		y += (x[i] - 100)*(x[i] - 100);
		l *= cos((x[i] - 100) / pow(i+1, 0.5));
	}
	y = y/4000-l+1;
	return y;
}
double H(double x, double t)
{
	if (x > t) return 1;
	else return 0;
}
double Gaussian()
{
	double a = double(rand()) / RAND_MAX, b = double(rand()) / RAND_MAX;
	return (cos(2 * pi*b))*pow(-2 * log(a), 0.5);
}
double garmma(double x)
{
	double temp[2], d = 1;
	temp[0] = 1 / (x + 1);
	unsigned i = 2;
	for (; d > pow(0.1, 8); i++)
	{
		temp[1] = temp[0] * i / (x + i);
		d = temp[0] * pow(i - 1, x) - temp[1] * pow(i, x);
		temp[0] = temp[1];
	}
	return temp[0] * pow(i, x);
}
class nest
{
public:
	vector<double> x;
	double quality,gar0,gar1;
	nest(){};
	nest(unsigned dim){ x.resize(dim, 0); gar0 = garmma(lamd); gar1 = garmma((1 + lamd) / 2); };
	void get_quality();
	void levy_fly();
	double similarity(nest &b);
};
void nest::get_quality()
{
	quality = 1 / (f(x) + 0.00001);
}
void nest::levy_fly()
{
	vector<double> y(x.size());
	double sigma = pow(gar0*sin(pi*lamd / 2) / (gar1*pow(2, (lamd - 1) / 2)), 1 / 2 * lamd), U, V, s, L;
	for (unsigned i = 0; i < x.size(); i++)
	{
		U = sigma*Gaussian();
		V = Gaussian();
		s = abs(U) / pow(abs(V), 1 / lamd);
		L = lamd*gar0*sin(pi*lamd / 2) / (pi*pow(s, 1 + lamd));
		y[i] = x[i] + alpha*L;
	}
	if (1 / (f(y) + 0.00001)>quality)
	{
		x = y;
		get_quality();
	}
}
class population
{
public:
	vector<nest> p;
	double e, b;
	population(double a, double c, unsigned n, unsigned dim){ b = a; e = c; p.resize(n, nest(dim)); init(0); get_quality(); };
	void init(unsigned h);
	void get_quality();
	void elim();
	void find();
};
void population::init(unsigned h)
{
	srand(h);
	for (unsigned i = 0; i < p.size(); i++)
	{
		for (unsigned j = 0; j < p[0].x.size(); j++)
		{
			p[i].x[j] = b + (e - b)*double(rand()) / RAND_MAX;
		}
	}
}
void population::get_quality()
{
	double m;
	unsigned k;
	for (unsigned i = 0; i < p.size(); i++)
	{
		p[i].get_quality();
		if (i == 0 || m<p[i].quality)
		{
			m = p[i].quality;
			k = i;
		}
	}
	nest temp = move(p[0]);
	p[0] = p[k];
	p[k] = temp;
}
void population::elim()
{
	vector<double> temp(p[0].x.size());
	for (unsigned i = 0; i < p.size(); i++)
	{
		int k[2];
		k[0] = rand() % p.size();
		k[1] = rand() % p.size();
		for (unsigned j = 0; j < p[0].x.size(); j++)
		{
			temp[j] = p[i].x[j] + af*(p[k[0]].x[j] - p[k[1]].x[j]);
		}
		p[i].get_quality();
		if (p[i].quality < 1 / (0.00001 + f(temp))) p[i].x = temp;
	}
}
void population::find()
{
	for (unsigned i = 0; i < p.size(); i++)
	{
		p[i].levy_fly();
	}
}
class calculator
{
public:
	unsigned n;
	double begin, end;
	population *p;
	nest best;
	calculator(population &x, unsigned in, double a, double b){ begin = a; end = b; n = in; p = &x; best = p->p[0]; };
	void advance();
	void iteration(ofstream &fout);
};
void calculator::advance()
{
	p->elim();
	p->find();
	p->get_quality();
}
void calculator::iteration(ofstream &fout)
{
	for (unsigned i = 0; i < n; i++)
	{
		if (i % 10000 == 0) srand(rand() % unsigned(clock()));
		/*if (i == 0)
		{
			for (unsigned j = 0; j < p->p[i].x.size();j++)cout << p->p[i].x[j] << '\t';
			cout << endl;
		}*/
		advance();
		if (i == 0 || best.quality < 1 / (f(p->p[0].x) + 0.00001))
		{
			best = p->p[0];
			best.get_quality();
		}
	}
	for (unsigned i = 0; i < best.x.size(); i++)
	{
		cout << "x" << i + 1 << "=" << best.x[i] << '\t';
	}
	cout << endl << "y=" << f(best.x) << endl;
	for (unsigned i = 0; i < best.x.size(); i++)
	{
		fout << best.x[i] << '\t';
	}
	fout << endl;
}
int main()
{
	clock_t start_time = clock();
	srand(unsigned(clock()));
	ofstream fout("result", ios::out | ios::trunc);
	fout.precision(12);
	unsigned pop_num = 200, it_num = 10000, dim = 4;
	double a = -600, b = 600;
	srand(unsigned(clock()));
	population p(a, b, pop_num, dim);
	calculator C(p, it_num, a, b);
	//for (unsigned i = 0; i < 100; i++)
	//{
		C.p->init(rand());
		C.iteration(fout);
	//}	
	fout.close();
	clock_t end_time = clock();
	cout << "time=" << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	getchar();
	return 0;
}
