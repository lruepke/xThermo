#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;
#include "H2ONaCl.h"
using namespace xThermal;
struct PROPS
{
    mxArray* T;
    mxArray* P;
    mxArray* X;
    mxArray* H;
    mxArray* S_l;
    mxArray* S_v;
    mxArray* S_h;
    mxArray* X_l;
    mxArray* X_v;
    mxArray* Rho_l;
    mxArray* Rho_v;
    mxArray* Rho_h;
    mxArray* H_l;
    mxArray* H_v;
    mxArray* H_h;
    mxArray* Rho;
    mxArray* phase;
};
void mGetMatrix(const mxArray *prhs, double **out, const char *varname, mwSize *out_m, mwSize *out_n)
{
    mwSize m, n;
    double *temp;

    m = mxGetM(prhs);
    n = mxGetN(prhs);

    if(!mxIsDouble(prhs)){
        char buff[256];
        sprintf(buff, "'%s' must be of type 'double'.\n", varname);
        mexErrMsgTxt(buff);
    }

    temp   = (double*)mxGetData(prhs);
    *out   = temp;
    *out_m = m;
    *out_n = n;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    string backend_name = "IAPS84";
    // check arguments
    if(nrhs!=3 && nrhs!=4){
        mexErrMsgTxt("Usage: [Rho_L, Rho_V] = X_VL(T, P, X)\n[Rho_L, Rho_V] = P_X_Critical(T, P, X, backend_name);");
    }
    if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) )
    {
        mexErrMsgTxt("The input T,P,X must be of type 'double'.");
    }
    if(nrhs==4)
    {
        int ind_backendName = nrhs-1;
        if (!mxIsChar(prhs[ind_backendName]))
        {
            mexErrMsgTxt("The second argument is used to specify which backend of H2O EOS will be used, it must be char, e.g. 'IAPS84', 'IAPWS95'. The default name is 'IAPS84'. ");
        }
        backend_name = mxArrayToString(prhs[ind_backendName]);
    }
    // mwSize nT = mxGetNumberOfElements(prhs[0]);
    double *TMat, *PMat, *XMat;
    mwSize m,n;
    mGetMatrix(prhs[0], &TMat, "T", &m, &n);
    mGetMatrix(prhs[1], &PMat, "P", &m, &n);
    mGetMatrix(prhs[2], &XMat, "X", &m, &n);
    double *prop_L_out, *prop_V_out;
    // plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    // prop_L_out = (double*)mxGetData(plhs[0]);
    // plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    // prop_V_out = (double*)mxGetData(plhs[1]);
    //calculate
    H2ONaCl::cH2ONaCl sw(backend_name);
    // cout<<"Use "<<backend_name<<" EOS of H2O"<<endl;

    const char* field_names[] = {"T", "P", "X", "H", "S_l", "S_v", "S_h", "X_l", "X_v", "Rho_l", "Rho_v", "Rho_h", "H_l", "H_v", "H_h", "Rho", "phase"};
    mwSize dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, dims, 17, field_names);
    // int name_T= mxGetFieldNumber(plhs[0], "T");
    mxArray* T_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* P_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* X_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* H_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* S_l_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* S_v_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* S_h_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* X_l_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* X_v_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* Rho_l_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* Rho_v_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* Rho_h_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* H_l_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* H_v_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* H_h_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* Rho_value = mxCreateDoubleMatrix(m, n, mxREAL);
    mxArray* phase_value = mxCreateDoubleMatrix(m, n, mxREAL);

    double *ptr_T = (double*)mxGetData(T_value);
    double *ptr_P = (double*)mxGetData(P_value);
    double *ptr_X = (double*)mxGetData(X_value);
    double *ptr_H = (double*)mxGetData(H_value);
    double *ptr_S_l = (double*)mxGetData(S_l_value);
    double *ptr_S_v = (double*)mxGetData(S_v_value);
    double *ptr_S_h = (double*)mxGetData(S_h_value);
    double *ptr_X_l = (double*)mxGetData(X_l_value);
    double *ptr_X_v = (double*)mxGetData(X_v_value);
    double *ptr_Rho_l = (double*)mxGetData(Rho_l_value);
    double *ptr_Rho_v = (double*)mxGetData(Rho_v_value);
    double *ptr_Rho_h = (double*)mxGetData(Rho_h_value);
    double *ptr_H_l = (double*)mxGetData(H_l_value);
    double *ptr_H_v = (double*)mxGetData(H_v_value);
    double *ptr_H_h = (double*)mxGetData(H_h_value);
    double *ptr_Rho = (double*)mxGetData(Rho_value);
    double *ptr_phase = (double*)mxGetData(phase_value);
    H2ONaCl::ThermodynamicProperty prop_cpp;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // mexPrintf("\ni: %" FMT_SIZE_T "u: %f \n", j+n*i, TMat[j + i*n]);
            mwSize ind = j + i*n;
            sw.prop_VL(TMat[ind], PMat[ind], XMat[ind], prop_cpp);
            // cout<<"T: "<<TMat[ind]<<"p: "<<PMat[ind]<<" Rhol: "<<prop_cpp.Rho_l<<" Xl: "<<prop_cpp.X_l<<endl;
            ptr_T[ind] = prop_cpp.state.T;
            ptr_P[ind] = prop_cpp.state.p;
            ptr_X[ind] = prop_cpp.state.X;
            ptr_H[ind] = prop_cpp.state.H;
            ptr_phase[ind] = prop_cpp.state.phase;
            ptr_S_l[ind] = prop_cpp.S_l;
            ptr_S_v[ind] = prop_cpp.S_v;
            ptr_S_h[ind] = prop_cpp.S_h;
            ptr_X_l[ind] = prop_cpp.X_l;
            ptr_X_v[ind] = prop_cpp.X_v;
            ptr_Rho_l[ind] = prop_cpp.Rho_l;
            ptr_Rho_v[ind] = prop_cpp.Rho_v;
            ptr_Rho_h[ind] = prop_cpp.Rho_h;
            ptr_H_l[ind] = prop_cpp.H_l;
            ptr_H_v[ind] = prop_cpp.H_v;
            ptr_H_h[ind] = prop_cpp.H_h;
            ptr_Rho[ind] = prop_cpp.Rho;
        }
    }
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "T"), T_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "P"), P_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "X"), X_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "H"), H_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "phase"), phase_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "S_l"), S_l_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "S_v"), S_v_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "S_h"), S_h_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "X_l"), X_l_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "X_v"), X_v_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "H_l"), H_l_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "H_v"), H_v_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "H_h"), H_h_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "Rho_l"), Rho_l_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "Rho_v"), Rho_v_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "Rho_h"), Rho_h_value);
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "Rho"), Rho_value);
}
