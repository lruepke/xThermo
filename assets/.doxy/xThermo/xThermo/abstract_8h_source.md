

# File abstract.h

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_templateClass**](dir_954facde36cda98a0027ddce7db85b3d.md) **>** [**abstract.h**](abstract_8h.md)

[Go to the documentation of this file](abstract_8h.md)


```C++

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
```


