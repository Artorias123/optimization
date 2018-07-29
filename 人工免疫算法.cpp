#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include <algorithm>
using namespace std;
#define gr 4294967296

double of(vector<double> x)
{
	double y=0;
	for (int i = 0; i < 2; i++)
	{
		y += (x[i] - 100)*(x[i] - 100);
	}
	y = y / 4000 - cos(x[0] - 100)* cos((x[1] - 100) / pow(2, 0.5)) + 1;
	return y;
}

class antibody
{
public:
	int n, antibody_concentration, repro_rate;
	vector<unsigned long> x;
	double antigen_affinity;
	antibody(int m=1);
	antibody(vector<unsigned long> y){ x = y; antibody_concentration = 0; };
	void input(vector<unsigned long> y){ x = y; antibody_concentration = 0; };
	~antibody(){};
	void init();
	void variation();
	void reproduction(antibody y, antibody z);
};
antibody::antibody(int m)
{ 
	n = m;
	antibody_concentration = 0;
	x.resize(n);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < 8; i++)
			x[j] += (rand() % 16)*pow(16, i);
	}
};
void antibody::init()
{
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < 7; i++)
			x[j] += (rand() % 16)*pow(16, i);
	}
}
void antibody::variation()
{
	for (int j = 0; j < n; j++)
	{
		unsigned i = rand() % 32;
		x[j] = x[j] ^ (1 << (i - 1));
	}
}
void antibody::reproduction(antibody y, antibody z)
{
	bool i;
	int k;
	for (int j = 0;j < n; j++)
	{
		unsigned i = rand() % 32;
		z.x[j] = z.x[j] ^ (1 << (i - 1));
		i = rand() % 2;
		k = rand() % 16+16;
		if (i) x[j] = z.x[j] | (y.x[j]>>k);
		else x[j] = z.x[j] | (y.x[j] << k);
	}
}
bool rcamp(antibody &a, antibody &b)
{
	return (a.repro_rate>b.repro_rate);
}

class population
{
public:
	int n;
	vector<antibody> indivial;
	population(int m,int dim=1);
	~population(){};
	void init();
	void get_antibody_concentration();
	void variation(double r);
	void select(double r,int h);
};
population::population(int m,int dim)
{
	n = m;
	indivial.resize(n, antibody(dim));
}
void population::init()
{
	for (int i = 0; i < n; i++)
	{
		indivial[i].init();
	}
}
void population::get_antibody_concentration()
{
	unsigned long p;
	for (int i = 0; i < n; i++)
	{
		indivial[i].antibody_concentration = 0;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			for (int k = 0; k < indivial[i].n; k++)
			{
				p = (indivial[i].x[k] ^ indivial[j].x[k]);
				for (int l = 0; l <10; l++)
				{
					if (!(p >> l))
					{
						indivial[i].antibody_concentration += 1;
						indivial[j].antibody_concentration += 1;
					}
				}
			}
		}
	}
}
void population::variation(double r)
{
	int d = ceil(n*r);
	for (int i = 0; i < d; i++)
	{
		indivial[rand() % n].variation();
	}
}
void population::select(double r,int h)
{
	int d = ceil(n*r), sr = 0;
	double s;
	vector<unsigned> re(d);
	vector<unsigned> rr(n,0);
	for (unsigned i = 0; i < n; i++)
	{
		sr += indivial[i].repro_rate;
		for (unsigned j = 0; j < i+1; j++)
		{
			rr[i] += indivial[i].repro_rate;
		}

	}
	for (int i = 0; i < d; i++)
	{
		re[i] = rand() % sr;
		for (unsigned j = 0; ; j++)
		{
			if (re[i]<rr[j])
			{
				re[i] = j;
				rand();
				break;
			}
			if (j==n-1)
			{
				re[i] = j;
				break;
			}
		}
	}
	sort(re.begin(), re.end(), less<unsigned>());
	for (int i = 0; i < d; i++)
	{
		indivial[i] = indivial[re[i]];
	}
}

class antigen
{
public:
	double a, l;
	antigen(double ia, double il){ a = ia; l = il; };
	~antigen(){};
	double f(vector<unsigned long> x);
	void match_antibody(population* p);
};
double antigen::f(vector<unsigned long> x)
{
	double y;
	vector<double> rx(x.size());
	for (unsigned i = 0; i < x.size(); i++)
	{
		rx[i] = a + l*x[i]/gr;
	}
	y = of(rx);
	return y;
}
void antigen::match_antibody(population* p)
{
	for (unsigned i = 0; i < p->n; i++)
	{
		p->indivial[i].antigen_affinity=1/(f(p->indivial[i].x)+0.01);
	}
}

class memory_library
{
public:
	vector<antibody> l;
	int n;
	int dim;
	memory_library(int m){ n = m; l.resize(n); };
	void get_seed(population* p);
};
bool acomp(antibody&a, antibody&b)
{
	return(a.antigen_affinity > b.antigen_affinity);
}
void memory_library::get_seed(population* p)
{
	sort(p->indivial.begin(), p->indivial.end(), acomp);
	for (int i = 0; p->indivial[i].antigen_affinity>l[n-1].antigen_affinity&&i<n; i++)
	{
		l[n-1] = p->indivial[i];
		sort(l.begin(), l.end(), acomp);
	}
	
}

class calculator
{
public:
	int iter;
	float vir_rate, pro_rate, diversity;
	population* p;
	antigen* a;
	memory_library* m;
	calculator(int it, float vr, float pr, float div,population& ip,memory_library& im,antigen& ia);
	void reproduction();
	void get_repro_rate();
	void get_anti_aff();
	void advance();
	void init();
	void iteration();
	antibody result();
};
calculator::calculator(int it, float vr, float pr, float div, population& ip, memory_library& im, antigen& ia)
{
	iter = it;
	vir_rate = vr;
	pro_rate = pr;
	diversity = div;
	p = &ip;
	m = &im;
	a = &ia;
}
void calculator::get_repro_rate()
{
	double sc = 0, sa = 0, b = 1 - diversity;
	for (unsigned i = 0; i < p->indivial.size(); i++)
	{
		sc += p->indivial[i].antibody_concentration;
		sa += p->indivial[i].antigen_affinity;
	}
	if (sc == 0) sc++;
	for (unsigned i = 0; i < p->indivial.size(); i++)
	{
		float a = diversity*(p->indivial[i].antigen_affinity / sa) + b*(p->indivial[i].antibody_concentration / sc);
		p->indivial[i].repro_rate = ceil(1000 * (diversity*(p->indivial[i].antigen_affinity / sa) + b*(p->indivial[i].antibody_concentration / sc)));
	}
}
void calculator::reproduction()
{
	int h = ceil(pro_rate*p->n), d =h +m->n;
	for (int i = h; i < d; i++)
	{
		p->indivial[i] = m->l[i-h];
	}
	for (int i = d; i <p->n; i++)
	{
		p->indivial[i].reproduction(p->indivial[rand() % d],p->indivial[rand() % d]);
	}
}
void calculator::advance()
{
	p->select(pro_rate,m->n);
	reproduction();
	p->variation(vir_rate);
	p->get_antibody_concentration();
	a->match_antibody(p);
	get_repro_rate();
	m->get_seed(p);
}
void calculator::init()
{
	srand(unsigned(clock()));
	p->init();
	p->get_antibody_concentration();
	a->match_antibody(p);
	get_repro_rate();
	m->get_seed(p);
	sort(m->l.begin(), m->l.end(), acomp);
}
void calculator::iteration()
{
	init();
	for (int i = 0; i < iter; i++)
	{
		advance();
	}
}
antibody calculator::result()
{
	return m->l[0];
}
class result
{
public:
	antibody ra;
	double a, l;
	vector<double> r;
	result(double ia, double il){ a = ia; l = il; };
	~result(){};
	void trans();
	void prt();
};
void result::trans()
{
	r.resize(ra.x.size());
	for (unsigned i = 0; i < r.size(); i++)
	{
		r[i] = a + l*ra.x[i] / gr;
	}
}
void result::prt()
{
	for (unsigned i = 0; i < r.size(); i++)
	{
		cout << r[i] << '\t';
	}
	cout << endl;
}
int main()
{
	clock_t start_time = clock();
	int pn=100, mm=10,it=10000,dim=2;
	float vr=0.6, pr=0.4, di=0.95,a=-600,h=1200;
	population p(pn, dim);
	memory_library m(mm);
	antigen f(a,h);
	calculator H(it,vr,pr,di,p,m,f);
	result r(a,h);
	H.iteration();
	r.ra = H.result();
	r.trans();
	double y = of(r.r);
	r.prt();
	cout << y << endl;
	clock_t end_time = clock();
	cout << "time=" << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	system("pause");
	return 0;
}