%module test_swig
%{
    #include "test.H"
%}

// %include std_string.i
// %include std_vector.i
// // ===== This is required for python API========
// namespace std {
//    %template(IntVector) vector<int>;
//    %template(DoubleVector) vector<double>;
//    %template(StringVector) vector<string>;
//    %template(ConstCharVector) vector<const char*>;
//    %template(UnsignedLongVector) vector<unsigned long int>;
// }
// // ========================================================================
// %inline %{ 
//     using namespace std;
// %}
// %include "typemaps.i"

namespace TEST_SWIG
{
    void test_print();
}