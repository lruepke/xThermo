#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;
#include "H2ONaCl.h"
using namespace xThermal;

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
    if(nrhs!=2 && nrhs!=3){
        mexErrMsgTxt("Usage: [X_L, X_V] = X_VL(T, P)\n[X_L, X_V] = P_X_Critical(T, P, backend_name);");
    }
    if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[0]))
    {
        mexErrMsgTxt("The T,P must be of type 'double'.");
    }
    if(nrhs==3)
    {
        int ind_backendName = nrhs-1;
        if (!mxIsChar(prhs[ind_backendName]))
        {
            mexErrMsgTxt("The second argument is used to specify which backend of H2O EOS will be used, it must be char, e.g. 'IAPS84', 'IAPWS95'. The default name is 'IAPS84'. ");
        }
        backend_name = mxArrayToString(prhs[ind_backendName]);
    }
    // mwSize nT = mxGetNumberOfElements(prhs[0]);
    double *TMat, *PMat;
    mwSize m,n;
    mGetMatrix(prhs[0], &TMat, "T", &m, &n);
    mGetMatrix(prhs[1], &PMat, "P", &m, &n);
    double *X_L_out, *X_V_out;
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    X_L_out = (double*)mxGetData(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    X_V_out = (double*)mxGetData(plhs[1]);
    //calculate
    H2ONaCl::cH2ONaCl sw(backend_name);
    // cout<<"Use "<<backend_name<<" EOS of H2O"<<endl;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            // mexPrintf("\ni: %" FMT_SIZE_T "u: %f \n", j+n*i, TMat[j + i*n]);
            mwSize ind = j + i*n;
            // cout<<"["<<i<<", "<<j<<"]"<<ind<<": "<<TMat[ind]<<endl;
            sw.X_VL(TMat[ind], PMat[ind], X_L_out[ind], X_V_out[ind]);
        }
    }
}
