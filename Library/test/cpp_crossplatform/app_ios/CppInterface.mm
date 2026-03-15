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
