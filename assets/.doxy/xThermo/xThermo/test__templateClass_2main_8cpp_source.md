

# File main.cpp

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_templateClass**](dir_954facde36cda98a0027ddce7db85b3d.md) **>** [**main.cpp**](test__templateClass_2main_8cpp.md)

[Go to the documentation of this file](test__templateClass_2main_8cpp.md)


```C++

#include "abstract.h"

#include "implement_A.h"
#include "implement_B.h"
#include "implement_AB.h"

#include <iostream>
using namespace std;


int main()
{
    cout<<"测试C++的泛型编程特性"<<endl;
    // IMPLEMENT_A::implement_A A;
    // cout<<A.calc_Tmax()<<endl;
    // IMPLEMENT_B::implement_B B;
    // cout<<B.calc_Tmax()<<endl;

    IMPLEMENT_AB::implement_AB A("A");
    cout<<A.calc_Tmax()<<endl;

    IMPLEMENT_AB::implement_AB B("B");
    cout<<B.calc_Tmax()<<endl;

    IMPLEMENT_AB::implement_AB C("C");
    cout<<C.calc_Tmax()<<endl;

    return 0;
}
```


