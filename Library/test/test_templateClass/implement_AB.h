#ifndef IMPLEMENT_AB_H
#define IMPLEMENT_AB_H

#include "abstract.h"
#include "implement_A.h"
#include "implement_B.h"
#include "Exception.h"

namespace IMPLEMENT_AB
{
    using namespace ABSTRACT;

    class implement_AB : public abstract
    {
    private:
        /* data */
    public:
        implement_AB(std::string which);
        ~implement_AB();
    private:
        abstract* m_which;
    public:
        double calc_Tmax();
    };

    double implement_AB::calc_Tmax()
    {
        return m_which->calc_Tmax();
    }
    
    implement_AB::implement_AB(std::string which)
    {
        if(which=="A")
        {
            STATUS("使用A实现");
            m_which = new IMPLEMENT_A::implement_A;
        }else if(which=="B")
        {
            STATUS("使用B实现");
            m_which = new IMPLEMENT_B::implement_B;
        }else
        {
            throw NotImplementedError("支持的选项为A,B其中之一"); 
        }
        
    }
    
    implement_AB::~implement_AB()
    {
        if(m_which)delete m_which;
    }
    
};

#endif