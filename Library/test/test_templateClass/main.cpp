
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