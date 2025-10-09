#include "demo.h"

double sum(double a, double b )
{
    return a + b;
}

void sum(std::vector<double> a, std::vector<double> b, std::vector<double>& result)
{
    result.clear();
    result.resize(a.size());
    for(size_t i = 0; i<a.size(); i++)
    {
        result[i] = a[i] + b[i];
    }
}