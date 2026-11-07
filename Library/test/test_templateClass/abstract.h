/**
 * @file abstract.h 
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief 定义顶层的抽象类
 * @version 0.1
 * @date 2022-03-26
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Exception.h"

#ifndef ABSTRACT_H
#define ABSTRACT_H

namespace ABSTRACT
{
    using namespace xThermal;

    class abstract
    {
    private:
        /* data */
    public:
        abstract(/* args */);
        ~abstract();
    public:
        virtual double calc_Tmax(){ throw NotImplementedError("calc_Tmax is not implemented for this backend"); };
        // virtual double calc_Tmin() = 0;
    };
    
    abstract::abstract(/* args */)
    {
    }
    
    abstract::~abstract()
    {
    }
    
};

#endif