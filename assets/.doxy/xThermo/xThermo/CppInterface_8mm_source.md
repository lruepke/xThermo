

# File CppInterface.mm

[**File List**](files.md) **>** [**app\_ios**](dir_e936cb0e42e32aaf0ec84f84716b9a91.md) **>** [**CppInterface.mm**](CppInterface_8mm.md)

[Go to the documentation of this file](CppInterface_8mm.md)


```C++
#import "CppInterface.h"
#include <example/example.h>
//#include "vector.H"

@interface CppInterface () 
{
    cExample* example;
}
@end

@implementation CppInterface

-(instancetype)init
{
    self = [super init];
    if (self) {
//        helloword();
        example = new cExample();
//        myFoo->PrintFoo();
        delete example;
        example = nullptr;
    }
    return self;
}

@end
```


