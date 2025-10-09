#ifndef IMPLEMENT_B_H
#define IMPLEMENT_B_H

#include "abstract.h"

namespace IMPLEMENT_B
{
    using namespace ABSTRACT;

    class implement_B : public abstract
    {
    private:
        /* data */
    public:
        implement_B(/* args */);
        ~implement_B();
    public:
        double calc_Tmax();
    };

    double implement_B::calc_Tmax()
    {
        return 50;
    }
    
    implement_B::implement_B(/* args */)
    {
    }
    
    implement_B::~implement_B()
    {
    }
    
};

#endif