#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>
using namespace std;
#include "H2ONaCl.h"
#include "IAPS84.h"
using namespace xThermal;
#define Kelvin 273.15

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

xThermal::PhaseRegion lookup_TorHPX(xThermal::cxThermal* m_pEOS, string m_fileName_LUT, const int& dim, const double& TorH, const double& P, const double& X, double* props, bool isCal, bool isCout=false) {
    xThermal::PhaseRegion phaseRegion;
    // check input parameters
    if (dim==2)
    {
        auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup();
        double xyz_min_target[dim];
        switch (pLUT->m_const_which_var) {
            case LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP:
            {
                // if(!CheckRange_TorH(pLUT->m_TorH,TorH, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of T is zero according to CONST_X_VAR_TorHP
                // if(!CheckRange_P(P, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of P is zero according to CONST_X_VAR_TorHP
                auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, TorH, P, isCal);
                phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                if (isCout)
                {
                    if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Input: T = "<<TorH-Kelvin<<" deg.C = "<<TorH<<" K, P = "<<P<<" Pa")
                    else STATUS("Input: H = "<<TorH/1000<<" kJ/kg = "<<TorH<<" J/kg, P = "<<P<<" Pa");
                    STATUS("Constant in the 2D LUT: X = "<<pLUT->m_constZ<<" kg/kg");
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP:
            {
                // if(!CheckRange_X(X, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of T is zero according to CONST_X_VAR_TorHP
                // if(!CheckRange_P(P, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of P is zero according to CONST_X_VAR_TorHP
                auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, X, P, isCal);
                phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                if (isCout)
                {
                    STATUS("Input: P = "<<P<<" Pa, X = "<<X<<" kg/kg");
                    if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Constant in the 2D LUT: T = "<<pLUT->m_constZ-Kelvin<<" deg.C = "<<pLUT->m_constZ<<" K")
                    else STATUS("Constant in the 2D LUT: T = "<<pLUT->m_constZ/1000<<" kJ/kg = "<<pLUT->m_constZ<<" J/kg");
                }
            }
                break;
            case LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH:
            {
                // if(!CheckRange_X(X, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of T is zero according to CONST_X_VAR_TorHP
                // if(!CheckRange_TorH(pLUT->m_TorH, TorH, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of P is zero according to CONST_X_VAR_TorHP
                auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, X, TorH);
                phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
                if (isCout)
                {
                    if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Input: T = "<<TorH-Kelvin<<" deg.C = "<<TorH<<" K, X = "<<X<<" kg/kg")
                    else STATUS("Input: H = "<<TorH/1000<<" kJ/kg = "<<TorH<<" J/kg, X = "<<X<<" kg/kg");
                    STATUS("Constant in the 2D LUT: P = "<<pLUT->m_constZ<<" Pa");
                }
            }
                break;
            default: ERROR("pLUT->m_const_which_var of 2D AMR-LUT is wrong: "<<pLUT->m_const_which_var<<", it is impossible. Check source code: ThermodynamicProperties cSWEOSarg::calculateSinglePoint_PTX(double P, double T_K, double X, bool isCout)");
                break;
        }
        if (isCout)
        {
            int i = 0;
            STATUS("Properties lookup result");
            cout<<"    "<<COLOR_BLUE<<"Phase region: "<<COLOR_PURPLE<<xThermal::map_phase2name[phaseRegion]<<COLOR_DEFAULT<<endl;
            for(auto &item : pLUT->m_map_props)
            {
                cout<<"    "<<COLOR_BLUE<<item.second.longName<<COLOR_PURPLE<<": "<<props[i]<<COLOR_DEFAULT<<" "<<item.second.unit<<endl;
                i++;
            }
        }
    }else if(dim==3)
    {
        auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup();
        // double props[pLUT->m_map_props.size()];
        double xyz_min_target[dim];
        // if(!CheckRange_TorH(pLUT->m_TorH, TorH, pLUT->m_xyz_min[0], pLUT->m_xyz_max[0], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of T is zero according to CONST_X_VAR_TorHP
        // if(!CheckRange_P(P, pLUT->m_xyz_min[1], pLUT->m_xyz_max[1], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of P is zero according to CONST_X_VAR_TorHP
        // if(!CheckRange_X(X, pLUT->m_xyz_min[2], pLUT->m_xyz_max[2], "in AMR-LUT "+m_fileName_LUT))exit(0);//index of P is zero according to CONST_X_VAR_TorHP
        auto *targetLeaf = m_pEOS->lookup(props, xyz_min_target, TorH, P, X);
        phaseRegion = targetLeaf->qData.leaf->user_data->phaseRegion_cell;
        if (isCout)
        {
            if(pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T)STATUS("Input: T = "<<TorH-Kelvin<<" deg.C = "<<TorH<<" K, P = "<<P<<" Pa, X = "<<X<<" kg/kg")
            else STATUS("Input: T = "<<TorH/1000<<" kJ/kg = "<<TorH<<" J/kg, P = "<<P<<" Pa, X = "<<X<<" kg/kg");
            int i = 0;
            STATUS("Properties lookup result");
            cout<<"    "<<COLOR_BLUE<<"Phase region: "<<COLOR_PURPLE<<xThermal::map_phase2name[phaseRegion]<<COLOR_DEFAULT<<endl;
            for(auto &item : pLUT->m_map_props)
            {
                cout<<"    "<<COLOR_BLUE<<item.second.longName<<COLOR_PURPLE<<": "<<props[i]<<COLOR_DEFAULT<<" "<<item.second.unit<<endl;
                i++;
            }
        }
    }
    return phaseRegion;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Input: (lutfilename, fluidName, T, P, [X], [backend of H2O])
    string backend_name = "IAPS84";
    bool isH2ONaCl = false;
    string fluidName = "H2ONaCl";
    string file_LUT;
    // check arguments
    if(nrhs!=4 && nrhs!=5 && nrhs!=6){
        mexErrMsgTxt("In put parameters doesn't meet the requirement of the function definition.\nUsage: props = lookup(lutfile, T, P, X, fluidName)\nprops = lookup(lutfile, T, P, X, fluidName, backend_name);\nprops = lookup(lutfile, T, P, fluidName);\nprops = lookup(lutfile, T, P, fluidName, backend_name);");
    }
    // 1. LUT file name
    if (!mxIsChar(prhs[0]))
    {
        mexErrMsgTxt("The first argument is used to specify filename of AMR-LUT, it must be char, e.g. 'lut_constT_XP_10.bin'. ");
    }else
    {
        file_LUT = mxArrayToString(prhs[0]);
    }
    //2. fluid name
    if (!mxIsChar(prhs[1]))
    {
        mexErrMsgTxt("The second argument is used to specify fluid name, it must be char, e.g. 'H2O-NaCl', 'H2O'. ");
    }else
    {
        fluidName = mxArrayToString(prhs[1]);
    }
    // 3. Input variables for the specific fluid
    if (fluidName=="H2O-NaCl")
    {
        if (nrhs==5) //input arguments are: lutfile, fluidname, T, P, X
        {
            if( !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]))
            {
                mexErrMsgTxt("The input T,P,X must be of type 'double'.");
            }
        }else if(nrhs==6) //input arguments are: lutfile, fluidname, T, P, X, backend of H2O
        {
            if( !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4]))
            {
                mexErrMsgTxt("The input T,P,X must be of type 'double'.");
            }
            int ind_backendName = nrhs-1;
            if (!mxIsChar(prhs[ind_backendName]))
            {
                mexErrMsgTxt("The last argument is used to specify which backend of H2O EOS will be used, it must be char, e.g. 'IAPS84', 'IAPWS95'. The default name is 'IAPS84'. ");
            }
            backend_name = mxArrayToString(prhs[ind_backendName]);
        }else
        {
            mexErrMsgTxt("The input arguments must be 5 or 6 for H2O-NaCl fluid, in order: [lut file name, fluid name, T, P, X, [Optional: backend of H2O] ]");
        }
    }else
    {
        if (nrhs==4 ) //input arguments must be: lutfile, fluidname, T, P for single component fluid, e.g., H2O
        {
            if( !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
            {
                mexErrMsgTxt("The input T,P must be of type 'double'.");
            }
        }else if(nrhs==5)
        {
            if( !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
            {
                mexErrMsgTxt("The input T,P must be of type 'double'.");
            }
            int ind_backendName = nrhs-1;
            if (!mxIsChar(prhs[ind_backendName]))
            {
                mexErrMsgTxt("The last argument is used to specify which backend of H2O EOS will be used, it must be char, e.g. 'IAPS84', 'IAPWS95'. The default name is 'IAPS84'. ");
            }
            backend_name = mxArrayToString(prhs[ind_backendName]);
        }
        else
        {
            mexErrMsgTxt("The input arguments must be 4 or 5, in order: [lut file name, fluid name, T, P, [Optional: backend of H2O] ]");
        }

    }

    // mwSize nT = mxGetNumberOfElements(prhs[0]);
    double *TMat, *PMat, *XMat;
    mwSize m,n;
    bool isCalculate = false; // whether call EOS to calculate properties on the need-refined cells. default is false.
    //1. load LUT
    LOOKUPTABLE_FOREST::EOS_ENERGY m_TorH = LOOKUPTABLE_FOREST::EOS_ENERGY_T;
    bool printStatus_lut = false;
    // xThermal::cxThermal* m_pEOS = new xThermal::H2ONaCl::cH2ONaCl(backend_name);
    xThermal::cxThermal* m_pEOS;
    if (isH2ONaCl)  // salt water: H2O-NaCl
    {
        m_pEOS = new xThermal::H2ONaCl::cH2ONaCl(backend_name);
        mGetMatrix(prhs[2], &TMat, "T", &m, &n);
        mGetMatrix(prhs[3], &PMat, "P", &m, &n);
        mGetMatrix(prhs[4], &XMat, "X", &m, &n);
    } else    // non-saltwater fluid
    {
        m_pEOS = new xThermal::PROST::cIAPS84();
        mGetMatrix(prhs[2], &TMat, "T", &m, &n);
        mGetMatrix(prhs[3], &PMat, "P", &m, &n);
        mGetMatrix(prhs[2], &XMat, "X", &m, &n); //dummy X as a placeholder for non-saltwater fluid
    }
    m_pEOS->loadLUT(file_LUT,printStatus_lut);
    const int dim = m_pEOS->get_dim_lut_lookup();
    // ====================== check if the input fluid name is the same as fluid name in the LUT file ---------
    std::string fluidName_in_LUT = m_pEOS->name();
    // if (fluidName != fluidName_in_LUT)
    // {
    //     if (!(fluidName=="H2O" && (fluidName_in_LUT == Name_Backend_IAPS84 || fluidName_in_LUT == Name_Backend_IAPWS95 || fluidName_in_LUT == Name_Backend_IAPWS95_CoolProp))) //special process for H2O
    //     {
    //         std::string errorinfo = "The input fluid name must be the same as fluid name in the AMR-LUT file.\nThe fluid name in the AMR-LUT file is '" + fluidName_in_LUT + "', but the input fluid name is '" + fluidName + "'\nPlease make sure the request thermodynamic properties are consistent with the input AMR-LUT !!!";
    //         mexErrMsgTxt(errorinfo.c_str());
    //     }
    // }
    //=========================================================================================================
    int num_props_per_node = 0;
    std::string name_props;
    std::vector<std::string> shortNames_prop_vector;
    std::vector<xThermal::propInfo> propsInfo;
    switch (dim) {
        case 2:
        {
            auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)m_pEOS->get_pLUT_lookup();
            num_props_per_node = (int)pLUT->m_map_props.size();
            std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
            // std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
            // if (pLUT->m_TorH != m_TorH) ERROR("The input AMR-LUT: "<<file_LUT<<" is a 2D table in "<<name_space<<" space, but your input TPX/HPX space is not match with that of the LUT file.");
            for(auto &item : pLUT->m_map_props)
            {
                propsInfo.push_back(item.second);
                name_props += "," + std::string(item.second.shortName)+ "(" + std::string(item.second.unit) + ")";
                shortNames_prop_vector.push_back(item.second.shortName);
            }
        }
            break;
        case 3:
        {
            auto* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)m_pEOS->get_pLUT_lookup();
            num_props_per_node = (int)pLUT->m_map_props.size();
            std::string name_space = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "TPX" : "HPX");
            std::string name_TorH = (pLUT->m_TorH == LOOKUPTABLE_FOREST::EOS_ENERGY_T ? "T" : "H");
            // if (pLUT->m_TorH != m_TorH) ERROR("The input AMR-LUT: "<<file_LUT<<" is a 3D table in "<<name_space<<" space, but your input TPX/HPX space is not match with that of the LUT file.");
            for(auto &item : pLUT->m_map_props)
            {
                propsInfo.push_back(item.second);
                name_props += "," + std::string(item.second.shortName)+ "(" + std::string(item.second.unit) + ")";
                shortNames_prop_vector.push_back(item.second.shortName);
            }
        }
            break;
        default: ERROR("Something is wrong in the AMR-LUT, because the dim is neither 2 nor 3, please check the LUT file: "<<file_LUT);
            break;
    }

    //
    xThermal::PhaseRegion phaseRegion;
    double props_per_node[num_props_per_node];
    propsInfo.push_back({"Phase","Phase region","-"});
    size_t N = m*n; //size of input T | P | X

    //========================= create matlab output matrix ==================
    const int n_props_per_node = shortNames_prop_vector.size();
    const char** prop_names_per_node = new const char* [n_props_per_node + 1]; //the last one is phase index
    std::vector<mxArray*> Props_value(n_props_per_node + 1); //the last one is phase index
    for (int i = 0; i < n_props_per_node; ++i) {
        prop_names_per_node[i] = shortNames_prop_vector[i].c_str();
        Props_value[i] = mxCreateDoubleMatrix(m, n, mxREAL);
    }
    Props_value[n_props_per_node] = mxCreateDoubleMatrix(m, n, mxREAL);
    prop_names_per_node[n_props_per_node] = "Phase";
    mwSize dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, dims, n_props_per_node + 1, prop_names_per_node);
    double* props[num_props_per_node + 1];//the last one is phase index
    for (int i = 0; i < num_props_per_node; ++i) {
        props[i]        = (double*)mxGetData(Props_value[mxGetFieldNumber(plhs[0], shortNames_prop_vector[i].c_str())]);
    }
    props[num_props_per_node]        = (double*)mxGetData(Props_value[mxGetFieldNumber(plhs[0], "Phase")]);
    //============================== LOOKUP ==========================================
    int j =0;
    for (int i = 0; i < N; ++i) {
        phaseRegion = lookup_TorHPX(m_pEOS, file_LUT, dim, TMat[i], PMat[i], XMat[i], props_per_node, isCalculate);
        for (int k = 0; k < num_props_per_node; ++k) {
            props[k][i] = props_per_node[k];
        }
        props[num_props_per_node][i] = phaseRegion;
    }

    // finish
    for (int i = 0; i < n_props_per_node; ++i) {
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], shortNames_prop_vector[i].c_str()), Props_value[i]);
    }
    mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "Phase"), Props_value[n_props_per_node]);
    delete[] prop_names_per_node;
    delete m_pEOS;
}
