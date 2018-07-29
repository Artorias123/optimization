#include<iostream>
#include<math.h>
#include<vector>
#include "stdlib.h"
#include <ctime>
using namespace std;
class nest
{
    public:
    vector<double> x;
    double quality,elim_rate;
    void random_init();
    void lavy_fly();
};
class population
{
    public:
    vector<nest> p;
    void init();
    void elim();
    void find();
};
class calculator
{
    public:
    unsigned ni;
};
