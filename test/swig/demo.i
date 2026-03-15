%module demo
%{
    #include "demo.h"
%}

%include std_string.i
%include std_vector.i
// ===== This is required for python API========
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(StringVector) vector<string>;
   %template(ConstCharVector) vector<const char*>;
   %template(UnsignedLongVector) vector<unsigned long int>;
}
// ========================================================================
%inline %{ 
    using namespace std;
%}
%include "typemaps.i"
%apply std::vector<double> *OUTPUT {std::vector<double>& result};

%include "demo.h"