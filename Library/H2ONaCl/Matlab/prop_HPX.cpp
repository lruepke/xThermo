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
    if(nrhs!=3 && nrhs!=4){
        mexErrMsgTxt("Usage: props = prop_HPX(H, P, X)\nprops = prop_HPX(H, P, X, backend_name);");
    }
    if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) )
    {
        mexErrMsgTxt("The input H,P,X must be of type 'double'.");
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
    double *HMat, *PMat, *XMat;
    mwSize m,n;
    mGetMatrix(prhs[0], &HMat, "H", &m, &n);
    mGetMatrix(prhs[1], &PMat, "P", &m, &n);
    mGetMatrix(prhs[2], &XMat, "X", &m, &n);
    double *prop_L_out, *prop_V_out;
    //calculate
    xThermal::cxThermal* pEOS = new xThermal::H2ONaCl::cH2ONaCl(backend_name);
    // H2ONaCl::cH2ONaCl sw(backend_name);
    // If add new property, please add names to field_names_vector and update the following stateArray
    std::vector<std::string> field_names_vector={"T", "P", "X", "H", "S_l", "S_v", "S_h", "X_l", "X_v", "Rho_l", "Rho_v", "Rho_h",
                                                 "H_l", "H_v", "H_h", "Mu_l", "Mu_v","Cp_l", "Cp_v", "Cp_h", "Rho", "Mu", "Cp", "phase"};
    //
    const int n_fields = field_names_vector.size();
    const char** field_names = new const char* [n_fields];
    std::vector<mxArray*> Fields_value(n_fields);
    for (int i = 0; i < n_fields; ++i) {
        field_names[i] = field_names_vector[i].c_str();
        Fields_value[i] = mxCreateDoubleMatrix(m, n, mxREAL);
    }
    mwSize dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, dims, n_fields, field_names);

    // call array version of the function
    size_t N = m*n;
    // !!!! The following names MUST be the same as field_names_vector
    xThermal::ThermodynamicPropertiesArray stateArray;
    stateArray.T        = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "T")]);
    stateArray.p        = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "P")]);
    stateArray.X        = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "X")]);
    stateArray.H        = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "H")]);
    stateArray.S_l      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "S_l")]);
    stateArray.S_v      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "S_v")]);
    stateArray.S_h      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "S_h")]);
    stateArray.X_l      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "X_l")]);
    stateArray.X_v      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "X_v")]);
    stateArray.Rho_l    = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Rho_l")]);
    stateArray.Rho_v    = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Rho_v")]);
    stateArray.Rho_h    = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Rho_h")]);
    stateArray.H_l      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "H_l")]);
    stateArray.H_v      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "H_v")]);
    stateArray.H_h      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "H_h")]);
    stateArray.Mu_l     = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Mu_l")]);
    stateArray.Mu_v     = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Mu_v")]);
    stateArray.Cp_l     = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Cp_l")]);
    stateArray.Cp_v     = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Cp_v")]);
    stateArray.Cp_h     = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Cp_h")]);
    stateArray.Rho      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Rho")]);
    stateArray.Mu       = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Mu")]);
    stateArray.Cp       = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Cp")]);
    stateArray.phase    = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "phase")]);
    pEOS->UpdateState_HPX(stateArray, N, HMat, PMat, XMat);

    for (int i = 0; i < n_fields; ++i) {
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], field_names_vector[i].c_str()), Fields_value[i]);
    }

    delete[] field_names;
}
