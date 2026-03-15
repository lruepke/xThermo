

# File test\_omp.cpp

[**File List**](files.md) **>** [**Library**](dir_5ad7f572bbca03234e8e621e192fc099.md) **>** [**test**](dir_fd7e605cf972a87e804527375e37bc17.md) **>** [**test\_omp.cpp**](test__omp_8cpp.md)

[Go to the documentation of this file](test__omp_8cpp.md)


```C++
#ifdef USE_OMP
#include "omp.h"
#endif
#include <iostream>
using namespace std;
void test_par()
{
#if USE_OMP == 1
// #pragma omp parallel
#pragma omp task
    cout<<omp_get_thread_num()<<endl;
#pragma omp task
    cout<<omp_get_thread_num()<<endl;
#pragma omp task
    cout<<omp_get_thread_num()<<endl;
#pragma omp task
    cout<<omp_get_thread_num()<<endl;
#pragma omp taskwait
    cout<<"wait"<<endl;
#endif
}
int main()
{
#if USE_OMP == 1
    omp_set_num_threads(8);
    int num = 0;
// #pragma omp parallel for
//     for (int i = 0; i < 10; ++i) {
//         cout<<i<<endl;
//         num = omp_get_num_threads();
//     }

    cout<<"num: "<<omp_get_num_threads()<<", "<<num<<endl;

#pragma omp parallel
    {
#pragma omp single
        {
            test_par();
        }

    };
#endif
    return 0;
}
```


