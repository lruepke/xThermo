

# File implement\_B.h

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_templateClass**](dir_954facde36cda98a0027ddce7db85b3d.md) **>** [**implement\_B.h**](implement__B_8h.md)

[Go to the documentation of this file](implement__B_8h.md)


```C++
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
```


