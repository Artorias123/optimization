#include "stdafx.h"
#include "iostream"
#include "math.h"
#include <vector>
#include <string>
using namespace std;
#define g_r (pow(5,0.5)-1)/2

double f(double x)//需要求极小值的函数
{
	double y;
	y = + 588600 / pow(10.8095222429746 - 4.21478541710781*cos(x - 2.09439333333333), 6)
		- 1079.1 / pow(10.8095222429746 - 4.21378541710781*cos(x - 2.09439333333333), 3)
		+ 481300 / pow(10.8095222429746 - 4.21378541710781*cos(x + 2.09439333333333), 6)
		- 1064.6 / pow(10.8095222429746 - 4.21378541710781*cos(x + 2.09439333333333), 3)
		+ 600800 / pow(10.8095222429746 - 4.21378541710781*cos(x), 6)
		- 1071.5 / pow(10.8095222429746 - 4.21378541710781*cos(x), 3);
	return y;
}

class Sig_peak//单峰区间
{
public:
	double begin, end;
	Sig_peak(double a, double b);
	~Sig_peak(){};
};
Sig_peak::Sig_peak(double a, double b)
{
	begin = a;
	end = b;
}

Sig_peak get_sigpeak(double a,double c, double h0)//前进-后退法搜索单峰区间
{
	double h = h0, b=a+h, d = h*0.1;
	bool k=0;
	for (int i = 0; b < c; i++)
	{
		if (f(a)>f(b))
		{
			k = 1;
		}
		else
		{
			if (f(a)>f(a+d)||k==1)
			{
				return Sig_peak(a-h, b);
			}
			else
			{
				k = 0;
			}
		}
		if (h<8*h0)
		{
			h *= 2;
		}
		 a = b;
		 b = b + h;
	}
	return Sig_peak(0, 0);
}
vector<Sig_peak> get_allpeak(double a, double b)//搜索待求区域[a,b]内全部单峰区间
{
	vector<Sig_peak> l;
	double h0, c;
	for (int j = 0; ; j++)
	{
		if (j==0) h0 = 0.05*(b - a);
		else (h0 = h0 / 5);
		c = a;
		vector<Sig_peak> temp;
		for (int i = 0;; i++)
		{
			temp.push_back(get_sigpeak(c, b, h0));
			if (temp[i].begin<temp[i].end)
			{
				c = temp[i].end;
			}
			else
			{
				temp.pop_back();
				break;
			}
		}
		if (l.size() != temp.size() || j == 0) l = temp;
		else break;
	}
	return l;
}

double get_min(double a, double d)//求单峰区间内极小值
{
	double c = g_r*(d - a) + a, b = a + g_r*(c - a);
	for (int i = 0; abs(c - b)>pow(10, -12); i++)
	{
		if (f(a)>f(b)&&f(b)>f(c))
		{
			a = b;
			b = c;
			c = d - g_r*(d - c);
		}
		else
		{
			d = c;
			c = b;
			b = a + g_r*(c - a);
		}
	}
	return (b + c) / 2;
}
vector<double> get_allmin(vector<Sig_peak> l)//求所有单峰区间内的极小值
{
	vector<double> m;
	int n = l.size();
	for (int i = 0; i < n; i++)
	{
		m.push_back(get_min(l[i].begin, l[i].end));
	}
	return m;
}

int _tmain(int argc, _TCHAR* argv[])
{
	vector<Sig_peak> l;
	vector<double> m;
	double y;
	int k=1;
	l = get_allpeak(0,10);
	m = get_allmin(l);
	int n = l.size();
	y = f(m[0]);
	for (int i = 0; i < n; i++)
	{
		cout << m[i] << endl;
		if (y>f(m[i]))
		{
			y = f(m[i]);
			k = i;
		}
	}
	cout << endl << y << endl << m[k] << endl;
	system("pause");
	return 0;
}