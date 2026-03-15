

# File implement\_A.h

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_templateClass**](dir_954facde36cda98a0027ddce7db85b3d.md) **>** [**implement\_A.h**](implement__A_8h.md)

[Go to the documentation of this file](implement__A_8h.md)


```C++
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
```


