#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;
#include "thermo.h"
#include "IAPS84.h"
#include "IAPWS95.h"
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
        mexErrMsgTxt("Usage: [props] = prop_water_TP(T, P)\n[props] = prop_water_TP(T, P, backend_name);");
    }
    if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) )
    {
        mexErrMsgTxt("The input T,P must be of type 'double'.");
    }
    if(nrhs==3)
    {
        int ind_backendName = nrhs-1;
        if (!mxIsChar(prhs[ind_backendName]))
        {
            mexErrMsgTxt("The second argument is used to specify which backend of H2O EOS will be used, it must be char, e.g. 'IAPS84', 'IAPWS95'. The default name is 'IAPS84'. ");
        }
        backend_name = mxArrayToString(prhs[ind_backendName]);
        if((backend_name!="IAPS84") && (backend_name != "IAPWS95"))
        {
            std::string info = "The input EOS name is not identified: "+backend_name+". The supported EOS name is one of  'IAPS84', 'IAPWS95'. Reset it to 'IAPS84'. ";
            mexWarnMsgTxt(info.c_str());
            backend_name = "IAPS84";
        }
    }
    // mwSize nT = mxGetNumberOfElements(prhs[0]);
    double *TMat, *PMat;
    mwSize m,n;
    mGetMatrix(prhs[0], &TMat, "T", &m, &n);
    mGetMatrix(prhs[1], &PMat, "P", &m, &n);

    //calculate
    xThermal::cxThermal* water;
    PROST::cIAPS84 iaps84;
    IAPWS95::cIAPWS95 iapws95;
    if (backend_name=="IAPS84")
    {
        water = &iaps84;
    }else if(backend_name=="IAPWS95")
    {
        water = &iapws95;
    } else
    {
        std::string info = "The input EOS name is not identified: "+backend_name;
        mexErrMsgTxt(info.c_str());
    }
    // If add new property, please add names to field_names_vector and update the following stateArray
    std::vector<std::string> field_names_vector={"T", "P", "H", "Rho", "Mu", "Cp", "phase"};
    //
    const int n_fields = field_names_vector.size();
    const char* field_names[n_fields];
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
    stateArray.H        = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "H")]);
    stateArray.Rho      = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Rho")]);
    stateArray.Mu       = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Mu")]);
    stateArray.Cp       = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "Cp")]);
    stateArray.phase    = (double*)mxGetData(Fields_value[mxGetFieldNumber(plhs[0], "phase")]);
    water->UpdateState_TPX(stateArray, N, TMat, PMat);

    for (int i = 0; i < n_fields; ++i) {
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], field_names_vector[i].c_str()), Fields_value[i]);
    }
}
