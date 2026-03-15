#include "H2ONaCl.h"
#include <fstream>
using namespace std;

void test_CriticalPX(int argc, char** argv);
void test_XL_VL(int argc, char** argv);
void help();
int main(int argc, char** argv)
{
    if(argc==1)help();
    int index_test = atoi(argv[1]);
    switch (index_test)
    {
    case 1:
        test_CriticalPX(argc, argv);
        break;
    case 2:
        test_XL_VL(argc, argv);
        break;
    default:
        help();
        break;
    }
}

void test_XL_VL(int argc, char** argv)
{
    H2ONaCl::cH2ONaCl sw("IAPWS95");
    // double T = 150 + 273.15, P=3.466540E5;
    double T =1.500000e+02 + 273.15, P =	3.492354E5;
    double XL=sw.XL_VL(T, P);
    printf("T = %.1f C, P = %.6f bar, XL = %.6e mol fraction\n", T - 273.15, P/1E5, sw.Wt2Mol(XL));
    printf("test VL\n");
}
void test_CriticalPX(int argc, char** argv)
{
    if(argc!=3)
    {
        printf("%sUsage%s: %s 1 [path of Electronic Annex EA-2 of Driesner(2007b)]\n",COLOR_GREEN,COLOR_DEFAULT,argv[0]);
        printf("%sUsage%s: %s 1 [Temperature value: %.2f-%.2f]\n",COLOR_GREEN,COLOR_DEFAULT,argv[0],H2ONaCl::T_MIN, H2ONaCl::T_MAX);
        exit(0);
    }
    H2ONaCl::cH2ONaCl sw("IAPWS95");
    double T = atof(argv[2]);
    string filename_EA2;
    if(T==0)
    {
        filename_EA2 = string(argv[2]);
        ifstream fin(filename_EA2);
        if(!fin)ERROR("Open file filed: "+filename_EA2);
        vector<double>T, P, X, P2, X2;
        double tmp,tmp_T, tmp_P, tmp_X,tmp_H, tmp_rho;
        string str_tmp;
        for(int i=0;i<6;i++)
        {
            getline(fin, str_tmp);
            // cout<<str_tmp<<endl;
        }
        while (!fin.eof())
        {
            fin>>tmp_T>>tmp_P>>tmp_X>>tmp_rho>>tmp_H;
            // cout<<tmp_T<<" "<<tmp_P<<" "<<tmp_X<<" "<<tmp_rho<<" "<<tmp_H<<endl; 
            T.push_back(tmp_T+273.15);
            P.push_back(tmp_P*1E5);
            X.push_back(tmp_X);   
        }
        fin.close();
        // calculate critical pressure and salinity
        sw.P_X_Critical(T, P2, X2);
        for(size_t i = 0; i<T.size();i++)
        {
            printf("T=%.2fK = %.2f deg.C, critical P err=%.6EPa, %.6EPa, critical X err=%f, %f wt%%NaCl\n", T[i],T[i]-273.15, P2[i], P[i], X2[i],X[i]);
        }
    }else
    {
        double P,X;
        sw.P_X_Critical(T, P, X);
        printf("T=%f K = %f deg.C, critical P=%.6e bar, critical X=%.6e wt%%NaCl = %.6e mol fraction\n", T, T-273.15, P/1E5, X, sw.Wt2Mol(X));
    }
}
void help()
{
    int test_index = 0;
    printf("%sPlease enter the following number to select test function.%s\n",COLOR_BLUE,COLOR_DEFAULT);
    test_index++; printf("%s%-3d%s: %sTest critical pressure and salinity%s\n",COLOR_GREEN, test_index,COLOR_DEFAULT,COLOR_PURPLE,COLOR_DEFAULT);
    test_index++; printf("%s%-3d%s: %sTest liquid composition at liquid branch of V+L coexistence%s\n",COLOR_GREEN, test_index,COLOR_DEFAULT,COLOR_PURPLE,COLOR_DEFAULT);
    exit(0);
}