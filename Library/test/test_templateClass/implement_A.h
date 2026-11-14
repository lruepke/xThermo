#ifndef IMPLEMENT_A_H
#define IMPLEMENT_A_H

#include "abstract.h"

namespace IMPLEMENT_A
{
    using namespace ABSTRACT;

    class implement_A : public abstract
    {
    private:
        /* data */
    public:
        implement_A(/* args */);
        ~implement_A();
    public:
        double calc_Tmax();
    };

    double implement_A::calc_Tmax()
    {
        return 20;
    }
    
    implement_A::implement_A(/* args */)
    {
    }
    
    implement_A::~implement_A()
    {
    }
    
};

#endif