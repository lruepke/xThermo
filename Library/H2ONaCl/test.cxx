/**
 * @file test.cxx
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Test H2ONaCl model
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "H2ONaCl.h"
using namespace std;
using namespace xThermal;

void test_constants(H2ONaCl::cH2ONaCl* sw);
void test_criticalCurve(H2ONaCl::cH2ONaCl* sw);
void test_HaliteLiquidus(H2ONaCl::cH2ONaCl* sw);
void debug_HaliteLiquidus_RhoL(H2ONaCl::cH2ONaCl* sw);
void debug_VL_props(H2ONaCl::cH2ONaCl* sw);
void test_phaseBoundaries2VTU(H2ONaCl::cH2ONaCl* sw);
void test_Props_TPX(H2ONaCl::cH2ONaCl* sw);
int lutGen_2D(int argc, char** argv);
int lutInfo(int argc, char** argv);
void test_T_VL(H2ONaCl::cH2ONaCl* sw);
void test_T_HPX(H2ONaCl::cH2ONaCl* sw);
void test_phaseRegion_HPX(xThermal::cxThermal* sw);
void test_viscosity(xThermal::cxThermal* sw);
void test_sliceP(H2ONaCl::cH2ONaCl* sw);
void test_sliceT(H2ONaCl::cH2ONaCl* sw);
void test_LUT_lookup(int argc, char** argv, xThermal::cxThermal& thermo);
void test_compressibility_multiphase(H2ONaCl::cH2ONaCl* sw, double S=0.5, double porosity=0.1, double rho_r=2000, double cp_r=1000);
int main(int argc, char** argv)
{
//    double T0=1000+273.15, P0=2.039018E8, X0=0.2;
//    STATUS("Test H2ONaCl");
   H2ONaCl::cH2ONaCl sw_84("IAPS84");
    H2ONaCl::cH2ONaCl sw_95("IAPWS95");
    H2ONaCl::cH2ONaCl sw_95_CoolProp("IAPWS95_CoolProp");
////    P0 = sw.P_VLH(T0); //(sw.P_VLH(T0) + sw.get_pWater()->Boiling_p(T0))/2.0;
//    double XV = sw.Wt2Mol(sw.XV_VL(T0, P0));
//    double rho = sw.Rho_phase(T0, 3.698853E5, sw.XL_VL(T0, 3.698853E5), Liquid);
//    double H, Cp;
//    sw.H_phase(T0, 3.698853E5, sw.XL_VL(T0, 3.698853E5), H, Cp, Liquid);
//     test_phaseRegion_HPX(&sw_84);
// test_T_HPX(&sw_84);
// test_T_VL(&sw);
// debug_HaliteLiquidus_RhoL(&sw);
// debug_VL_props(sw);
// test_phaseBoundaries2VTU(sw);
// test_Props_TPX(&sw_84);
//  lutGen_2D(argc, argv);
//    lutInfo(argc,argv);
//    test_viscosity(&sw);
//     test_sliceP(&sw_84);
//     test_sliceT(&sw_84);
// ThermodynamicPropertiesArray propsArray;
// propsArray.create(1);
// propsArray.T[0]=580.3 + 273.15;
// propsArray.p[0] = 955.6E5;
// propsArray.X[0] = 0.1;
// xThermal::cxThermal* eos = new xThermal::H2ONaCl::cH2ONaCl("IAPS84");
// eos->UpdateState_TPX(propsArray, 1, propsArray.T, propsArray.p, propsArray.X );
//     test_LUT_lookup(argc, argv, sw_84);

    // double p0 = 3.93188e+06, X0 = 0.3;
    // double h0 = 7.29702e+07/41.5534;
    // ThermodynamicProperties props;
    // sw_84.UpdateState_TPX(props, 538.431, p0, X0);
    // std::cout<<props<<endl;
    // sw_84.UpdateState_HPX(props, h0, p0, X0); //532.163
    // std::cout<<props<<endl;

    double p0 = 1.29884e+07, T0 = 631.905;
    double X0 = 0.349809;
    ThermodynamicProperties props;
    sw_84.UpdateState_TPX(props, T0, p0, X0);
    cout<<props<<endl;

    // double T0=438.76642135933298, p0=774575.32599110622, X0=0.42121320366209464;
    // ThermodynamicProperties props;
    // sw_84.UpdateState_TPX(props, T0, p0, X0);
    // cout<<props<<endl;
    // sw_84.UpdateState_TPX(props, T0, p0+1E5, X0);
    // cout<<props<<endl;

    // test_compressibility_multiphase(&sw_84);

    return 0;
}

void test_compressibility_multiphase(H2ONaCl::cH2ONaCl* sw, double S, double porosity, double rho_r, double cp_r)
{
    double T0 = 250 + 273.15, p0 = 39E5, X0 = 0.01;
    ThermodynamicProperties props;

    // 1. V+L
    X0 = sw->XL_VL(T0, p0);
    sw->UpdateState_TPX(props, T0, p0, X0);
    cout<<props<<endl;
    printf("V+L compressibility: %.2E per. bar\n", props.IsothermalCompressibility*1E5);
    printf("V+L compressibility(L): %.2E per. bar\n", props.IsothermalCompressibility_l*1E5);
    printf("V+L compressibility(V): %.2E per. bar\n", props.IsothermalCompressibility_v*1E5);

    // 2. l+H
    X0 = sw->X_HaliteLiquidus(T0, p0)*2;
    sw->UpdateState_TPX(props, T0, p0, X0);
    cout<<props<<endl;
    printf("L+H compressibility: %.2E per. bar\n", props.IsothermalCompressibility*1E5);
    printf("L+H compressibility(L): %.2E per. bar\n", props.IsothermalCompressibility_l*1E5);

    // 3. V+H
    T0=400 + 273.15;
    X0 = 0.3; //sw->X_VH(T0, p0);
    sw->UpdateState_TPX(props, T0, p0, X0);
    cout<<props<<endl;
    printf("V+H compressibility: %.2E per. bar\n", props.IsothermalCompressibility*1E5);
    printf("V+H compressibility(V): %.2E per. bar\n", props.IsothermalCompressibility_v*1E5);

    // 4. V+L+H
    ThermodynamicProperties props_vlh;
    p0 = sw->P_VLH(T0);
    double H = 2.14155e+06; //2.0E6;
    X0=0.6;
    sw->UpdateState_HPX(props_vlh, H, p0, X0);
    double Hmin, Hmax;
    sw->HminHmax_VLH_triangle(props_vlh.H_v, props_vlh.H_l, props_vlh.H_h, props_vlh.X_v, props_vlh.X_l, X0, Hmin, Hmax);
    double L = Hmax - Hmin; //latent heat
    double dPdT = sw->dPdT_VLH(props_vlh.T);
    double compressibility_vlh = ((1-porosity)*rho_r*cp_r + porosity*(props_vlh.S_l*props_vlh.Rho_l*props_vlh.Cp_l + props_vlh.S_v*props_vlh.Rho_v*props_vlh.Cp_v + props_vlh.S_h*props_vlh.Rho_h*props_vlh.Cp_h)) * (props_vlh.Rho_l - props_vlh.Rho_v) /(porosity * L * dPdT * props_vlh.Rho_l*props_vlh.Rho_v);
    cout<<props_vlh<<endl;
    printf("Hmin = %.4E, Hmax = %.4E, L = %f, dPdT = %.4E\n", Hmin, Hmax, L, dPdT);
    printf("Liquid phase compressibility: %.2E per. bar\n", props_vlh.IsothermalCompressibility_l*1E5);
    printf("Vapor phase compressibility: %.2E per. bar\n", props_vlh.IsothermalCompressibility_v*1E5);
    printf("Three phase compressibility(Grant): %.3f per. bar\n", compressibility_vlh*1E5);

    cout<<"P_vlh: "<<p0<<endl;
}

void test_LUT_lookup(int argc, char** argv, xThermal::cxThermal& thermo)
{
    using namespace std;
    // thermo.loadLUT("/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/hydrothermalfoam/tutorials/test/HydrothermalSinglePhaseDarcyFoam_xThermal_LUT/AMR_LUT_H2ONaCl/lut_TPX_9.bin");
    thermo.loadLUT("lut_constT_XP_10.bin");
   //simple lookup
   LOOKUPTABLE_FOREST::LookUpTableForest_2D* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)thermo.get_pLUT_lookup();
   double * props = new double[pLUT->m_map_props.size()];
   double xyz_min_target[3];
   double T0=600+273.15, p0=10E5, X0=0.3;
   // 1. lookup
   thermo.lookup(props, xyz_min_target,X0, p0, true);
   vector<string> propNames;
   for(auto & item : pLUT->m_map_props)propNames.push_back(item.second.shortName);
   cout<<"lut T "<<T0<<", p "<<p0<<", X "<<X0<<", rho: "<<props[0]<<endl;
    for (int i = 0; i < pLUT->m_map_props.size(); ++i) {
        cout<<"  "<<propNames[i]<<": "<<props[i]<<endl;
    }
   // 2. calculate
   ThermodynamicProperties prop;
   thermo.UpdateState_TPX(prop, T0, p0, X0);
   cout<<prop<<endl;

   exit(0);

    // generate data
    int N =60;
    double Tmin = 1+273.15,Tmax = 1000+273.15, pmin = 5.1E5, pmax=500E5;
    std::vector<double>T = thermo.linspace(Tmin, Tmax, N);
    std::vector<double>p = thermo.linspace(pmin, pmax, N);

    STATUS("Start looking up ... ");
    ThermodynamicProperties prop_lookup;
//    vector<H2ONaCl::PhaseRegion> pb_lookup(x.size(), H2ONaCl::SinglePhase_L);
    int ind = 0;
    time_t start = clock();
    int dim_lut = 2;
    switch (dim_lut)
    {
        case 2:
        {
            const int dim = 2;
            LOOKUPTABLE_FOREST::LookUpTableForest_2D* pLUT = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)thermo.get_pLUT_lookup();
            LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
            FILE* fp = NULL;
            fp = fopen("lookup_result.csv", "w");
            if(!fp)ERROR("Open file failed: lookup_result.csv");
            FILE* fp_eos = NULL;
            fp_eos = fopen("lookup_result_eos.csv", "w");
            if(!fp_eos)ERROR("Open file failed: lookup_result_eos.csv");
            double* props = new double[pLUT->m_map_props.size()];
            ThermodynamicProperties prop;
            double xyz_min_target[dim];
            fprintf(fp, "Phase region\tNeed_refine");
            int ind_Rho = -1, ind_Cp = -1, ind_Mu = -1, ind_kappa = -1, ind_beta = -1;
            int i = 0;
            for (auto &m : pLUT->m_map_props)
            {
                fprintf(fp, "\t%s", m.second.longName);
                switch (m.first) {
                    case Update_prop_Rho:
                        ind_Rho = i;
                        break;
                    case Update_prop_Cp:
                        ind_Cp = i;
                        break;
                    case Update_prop_Mu:
                        ind_Mu = i;
                        break;
                    case Update_prop_IsothermalCompressibility:
                        ind_kappa = i;
                        break;
                    case Update_prop_IsobaricExpansivity:
                        ind_beta = i;
                        break;
                }
                i++;
            }
            // check if all the required properties are contained in the LUT
            if(ind_Rho<0) ERROR("The Rho is not found in the input AMR-LUT");
            if(ind_Cp<0) ERROR("The Cp is not found in the input AMR-LUT");
            if(ind_Mu<0) ERROR("The Mu is not found in the input AMR-LUT");
            if(ind_kappa<0) ERROR("The kappa is not found in the input AMR-LUT");
            if(ind_beta<0) ERROR("The beta is not found in the input AMR-LUT");

            fprintf(fp, "\n");
            fprintf(fp_eos, "Bulk density\tBulk specific heat\tDynamic viscosity\tIsothermal compressibility\tIsobaric expansivity\n");
            for (size_t i = 0; i < T.size(); i++)
            {
                for (int jj = 0; jj < p.size(); ++jj) {
                    targetLeaf = thermo.lookup(props, xyz_min_target, T[i], p[jj], true);
                    fprintf(fp, "%d\t%d", targetLeaf->qData.leaf->user_data->phaseRegion_cell, targetLeaf->qData.leaf->user_data->need_refine);
                    for (size_t j = 0; j < pLUT->m_map_props.size(); j++)
                    {
                        fprintf(fp, "\t%.8E", props[j]);
                    }
                    fprintf(fp, "\n");
                    if(targetLeaf->qData.leaf->user_data->phaseRegion_cell == MixPhaseRegion)ind++;
                //    print eos result
                    thermo.UpdateState_TPX(prop, T[i],p[jj]);
                    fprintf(fp_eos, "%.8E/%.8E\t%.8E/%.8E\t%.8E/%.8E\t%.8E/%.8E\t%.8E/%.8E\n",
                            prop.Rho, props[ind_Rho], prop.Cp, props[ind_Cp], prop.Mu, props[ind_Mu],
                            prop.IsothermalCompressibility, props[ind_kappa], prop.IsobaricExpansivity, props[ind_beta]);
                }

            }
            fclose(fp);
            fclose(fp_eos);
            delete[] props;
        }
            break;
        default:
            break;
    }
    printf("All %d/%ld (%.2f %%) points close to phase boundary.\n", ind, T.size(), ind/(double)T.size()*100);
    STATUS_time("Searching done", clock() - start);
    // // write result to file
    // STATUS("Writting results to file");
    // ofstream fout("lookup_result.csv");
    // if(!fout)ERROR("Open file failed: lookup_result.csv");
    // fout<<"x\ty\tz\tphase region index\tphase region name"<<endl;
    // for (size_t i = 0; i < x.size(); i++)
    // {
    //     fout<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<"\t"<<pb_lookup[i]<<"\t"<<sw.m_phaseRegion_name[pb_lookup[i]]<<endl;
    // }

    // fout.close();
}

void test_sliceT(H2ONaCl::cH2ONaCl* sw)
{
    PhaseRegion_Slice slice= sw->Slice_constT(100);
    // std::vector<double> X = slice.regions["V+H"][0].X;
    // std::vector<double> T = slice.regions["V+H"][0].T;
    // for (int i = 0; i < 10; ++i) {
    //     cout<<i<<endl;
    //     // sw->Triangulation(X, T, sw->mean_vector(X), sw->mean_vector(T), 0.1, 1);
    //     for (int j = 0; j < X.size(); ++j) {
    //         X[j] = log10(slice.regions["V+H"][0].X[j]);
    //     }
    //     // X=slice.regions["V+L"][1].X;
    //     sw->Triangulation(X, T, sw->mean_vector(X), sw->mean_vector(T), 1, 1);
    // }
    std::cout<<slice<<endl;
}
void test_sliceP(H2ONaCl::cH2ONaCl* sw)
{
    PhaseRegion_Slice slice= sw->Slice_constP(250E5);
    // std::vector<double> X = slice.regions["V+H"][0].X;
    // std::vector<double> T = slice.regions["V+H"][0].T;
    // for (int i = 0; i < 10; ++i) {
    //     cout<<i<<endl;
    //     // sw->Triangulation(X, T, sw->mean_vector(X), sw->mean_vector(T), 0.1, 1);
    //     for (int j = 0; j < X.size(); ++j) {
    //         X[j] = log10(slice.regions["V+H"][0].X[j]);
    //     }
    //     // X=slice.regions["V+L"][1].X;
    //     sw->Triangulation(X, T, sw->mean_vector(X), sw->mean_vector(T), 1, 1);
    // }
    std::cout<<slice<<endl;
}

void test_viscosity(xThermal::cxThermal* sw)
{
    ThermodynamicProperties props;
    double T,P, X;
    // T= 900+273.15, P = 300E5, X=0.1;
    // T = 1 + 273.15, P = 200E5, X=0.1;
    T = 1 + 273.15; P = 100E5; X = 0.7;
    sw->UpdateState_TPX(props, T,P, X);
    double mu_l = props.Mu_l;
    double mu_v = props.Mu_v;
    double mu = props.Mu;
    PhaseRegion phase = props.phase;
    cout<<"phase: "<<sw->phase_name(phase)<<", mu: "<<mu<<", mu_l: "<<mu_l<<", mu_v: "<<mu_v<<endl;
}

void test_phaseRegion_HPX(xThermal::cxThermal* sw)
{
    ThermodynamicProperties props2= sw->UpdateState_TPX(373.15, 200E5, 0.2);
    cout<<props2.dRhodP<<" "<<props2.dRhodT<<endl;

    exit(0);

    double H = 4.5E6, X = 0.1, P = 600E5;
    ThermodynamicProperties props;
    sw->UpdateState_HPX(props,H,P,X);
    PhaseRegion phase = props.phase;
    cout<<sw->phase_name(phase)<<endl;
    std::vector<double> Hvec = {H};
    std::vector<double> Pvec = {P};
    std::vector<double> Xvec = {X};
    ThermodynamicPropertiesVector aa = sw->UpdateState_HPX(Hvec, Pvec, Xvec);
    std::cout<<aa<<endl;
    cout<<sw->phase_name((PhaseRegion)aa.phase[0])<<endl;
}

void test_T_HPX(H2ONaCl::cH2ONaCl* sw)
{
    double T = 1 + 273.15;
    double P = 10E5;
    double X = 1E-8;
    ThermodynamicProperties props;
    sw->UpdateState_TPX(props, T,P,X);
    cout<<"T = "<<T<<", p = "<<P<<", X="<<X<<", H="<<props.H<<endl;
    P = 22173339.843750;
    X = 0.000010;
    double T1, T2;
    sw->T_VLH_P0(P, T1, T2);
    double T_ = sw->T_HPX(1974023.437500, P, X);
    cout<<"T_ = "<<T_<<endl;
}

void test_T_VL(H2ONaCl::cH2ONaCl* sw)
{
    double T = H2ONaCl::T_MIN_VLH;
    // double P = 100E5;
    double p_vlh = sw->P_VLH(T);
    double p_boil = sw->get_pWater()->Boiling_p(T);
    double P = (p_vlh + p_boil)/2.0;
    double T_crit;
    sw->T_Critical(P, T_crit);
    double xl_vl = sw->XL_VL(T,P);
    double xv_vl = sw->XV_VL(T,P);
    double Hv_VL;
    sw->H_phase(T,P,xv_vl,Hv_VL, Vapor);
    double T_L = sw->T_VL_L(P, xl_vl, T_crit, sw->Tmax());
    double T_V = sw->T_VL_V(P, xv_vl, sw->Tmin(), sw->Tmax());
    cout<<"T = "<<T<<", P="<<P<<", XL = "<<xl_vl<<", xv_vl: "<<xv_vl<<", inverse T_L = "<<T_L<<", inverse T_V = "<<T_V<<", Hv: "<<Hv_VL<<endl;
    cout<<"xv_2 "<<sw->XV_VL(T_V, P)<<", p:bar: "<<P/1E5<<endl;
}

int lutGen_2D(int argc, char** argv)
{
    if(argc!=11)
    {
        STATUS("Usage: "+string(argv[0])+" [T|H] [T|H|P|X] [z0] [xmin] [xmax] [ymin] [ymax] [min_level] [max_level] [nThread]");
        STATUS("Unit: T[deg.C], P[bar], X[wt. NaCl, 0-1], H[MJ/kg]");
        STATUS_color("Example 1: T X 0.2 1 1000 5      600  4   7   8 (constX, var TP)", COLOR_BLUE);
        STATUS_color("Example 2: T P 200 0.001 1 1     1000 4   7   8 (constP, var XT)", COLOR_BLUE);
        STATUS_color("Example 3: T T 100 0.001 1 5     600  4   7   8 (constT, var XP)", COLOR_BLUE);
        STATUS_color("Example 5: H X 0.2 0.1 3.9 5     600  4   7   8 (constX, var HP)", COLOR_PURPLE);
        STATUS_color("Example 6: H P 200 0.001 1 0.1   3.9  4   7   8 (constP, var XH)", COLOR_PURPLE);
        STATUS_color("Example 7: H H 100 0.001 1 5     600  4   7   8 (constH, var XP)", COLOR_PURPLE);
        return 0;
    }
    string TorH = string(argv[1]);
    string const_var = string(argv[2]);
    double constZ = atof(argv[3]);
    double xmin = atof(argv[4]);
    double xmax = atof(argv[5]);
    double ymin = atof(argv[6]);
    double ymax = atof(argv[7]);
    int min_level = atoi(argv[8]);
    int max_level = atoi(argv[9]);
    int n_threads = atoi(argv[10]);
    H2ONaCl::cH2ONaCl eos("IAPWS95");
#if USE_OMP == 1
    eos.set_num_threads(n_threads);
#endif
    if(min_level > max_level)
    {
        WARNING("min_level > max_level: " + to_string(min_level)+" "+to_string(max_level));
        int tmp = min_level;
        min_level = max_level;
        max_level = tmp;
    }
    if(min_level<0)
    {
        WARNING("min_level<0: "+ to_string(min_level));
        min_level = 0;
    }
    if(max_level>(MAX_FOREST_LEVEL-3))
    {
        WARNING("min_level is too big: " + to_string(min_level));
        min_level = MAX_FOREST_LEVEL - 3;
    }
    if (TorH == "T")
    {
        if(const_var == "X")
        {
            eos.createLUT_2D(xmin + 273.15, xmax + 273.15, ymin * 1E5, ymax * 1E5, constZ, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_CommonlyUsedTPX);
            eos.save_lut_to_binary("lut_constX_TP_"+std::to_string(max_level)+".bin");
            eos.save_lut_to_vtk("lut_constX_TP_"+std::to_string(max_level)+".vtu");
        }else if(const_var == "P")
        {
            eos.createLUT_2D(xmin, xmax, ymin + 273.15, ymax+273.15, constZ * 1E5, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_CommonlyUsedTPX);
            eos.save_lut_to_binary("lut_constP_XT_"+std::to_string(max_level)+".bin");
            eos.save_lut_to_vtk("lut_constP_XT_"+std::to_string(max_level)+".vtu");
        }else if(const_var == "T")
        {
            eos.createLUT_2D(xmin, xmax, ymin*1E5, ymax*1E5, constZ + 273.15, LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_CommonlyUsedTPX);
            eos.save_lut_to_binary("lut_constT_XP_"+std::to_string(max_level)+".bin");
            eos.save_lut_to_vtk("lut_constT_XP_"+std::to_string(max_level)+".vtu");
        }
        else
        {
            ERROR("If the first arg is T, the second arg must be one of T, P, X");
        }
    }else if (TorH == "H")
    {
        if(const_var == "X")
        {
            eos.createLUT_2D(xmin*1E6, xmax*1E6, ymin*1E5, ymax*1E5, constZ, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_CommonlyUsedHPX);
            eos.save_lut_to_binary("lut_constX_HP_"+std::to_string(max_level)+".bin");
            eos.save_lut_to_vtk("lut_constX_HP_"+std::to_string(max_level)+".vtu");
        }else if(const_var == "P")
        {
            eos.createLUT_2D(xmin, xmax, ymin*1E6, ymax*1E6, constZ*1E5, LOOKUPTABLE_FOREST::CONST_P_VAR_XTorH, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_CommonlyUsedHPX);
            // WAIT("createLUT_2D");
            eos.save_lut_to_binary("lut_constP_XH_"+std::to_string(max_level)+".bin");
            eos.save_lut_to_vtk("lut_constP_XH_"+std::to_string(max_level)+".vtu");
        }else if(const_var == "H")
        {
            eos.createLUT_2D(xmin, xmax, ymin*1E5, ymax*1E5, constZ*1E6, LOOKUPTABLE_FOREST::CONST_TorH_VAR_XP, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_CommonlyUsedHPX);
            eos.save_lut_to_binary("lut_constH_XP_"+std::to_string(max_level)+".bin");
            eos.save_lut_to_vtk("lut_constH_XP_"+std::to_string(max_level)+".vtu");
        }
        else
        {
            ERROR("If the first arg is H, the second arg must be one of H, P, X");
        }
    }else
    {
        ERROR("The first arg must be T or H.");
    }

    return 0;
}

int lutInfo(int argc, char** argv)
{
    H2ONaCl::cH2ONaCl eos("IAPS84");
    eos.getLutInfo(argv[1]);
//    if(argc!=2)ERROR("Usage: lutinfo myLUT.bin");
//    // H2ONaCl::cH2ONaCl sw;
//    // sw.loadLUT(argv[1]);
//    int m_dim_lut = LOOKUPTABLE_FOREST::get_dim_from_binary(argv[1]);
//
//    switch (m_dim_lut)
//    {
//        case 2:
//        {
//            LOOKUPTABLE_FOREST::LookUpTableForest_2D* m_pLUT_lookup = (LOOKUPTABLE_FOREST::LookUpTableForest_2D*)(new LOOKUPTABLE_FOREST::LookUpTableForest_2D(argv[1], NULL));
//            printf("const %d, %f\n",m_pLUT_lookup->m_const_which_var, m_pLUT_lookup->m_constZ);
//        }
//            break;
//        case 3:
//        {
//            LOOKUPTABLE_FOREST::LookUpTableForest_3D* m_pLUT_lookup = (LOOKUPTABLE_FOREST::LookUpTableForest_3D*)(new LOOKUPTABLE_FOREST::LookUpTableForest_3D(argv[1], NULL));
//        }
//            break;
//        default:
//        ERROR("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+string(argv[1]));
//            break;
//    }
    return 0;
}

void test_Props_TPX(H2ONaCl::cH2ONaCl* sw)
{
    double T,P,X;
    xThermal::ThermodynamicProperties props;
    T = 1273.15, P = 250E5, X=1;
    sw->UpdateState_TPX(props, T,P,X);
    double H = props.H;

    int nT = 10;
    double Tmin = 2 + 273.15, Tmax = 650 + 273.15;
    double dT = (Tmax - Tmin)/(nT-1);
    P = 5.057798e+02 * 1E5;
    X = 0.1;
    time_t start = clock();

    for (int i = 0; i < nT; ++i) {
        sw->UpdateState_TPX(props, Tmin + dT*i, P, X);
        std::cout<<props.Rho<<" "<<props.H<<" "<<sw->phase_name(props.phase)<<std::endl;
    }
    time_t end = clock();
    STATUS_time("UpdateState_TPX",end-start);
}

void test_phaseBoundaries2VTU(H2ONaCl::cH2ONaCl* sw)
{
    PhaseBoundaries pb = sw->calc_PhaseBoundaries("linear",1, 0.1);
    sw->normalizePhaseBoundaries(pb);
    sw->writePhaseBoundaries2VTU("pb_linear", pb, 1, 1, 1);
}

void debug_HaliteLiquidus_RhoL(H2ONaCl::cH2ONaCl* sw)
{
    double  T=0;
    //1.  T=3.000000e+01, Rhol=nan, see also mmc3 of Driesner2007b
     T=3.000000e+01 + 273.15;

    //2.  High T low P: Hl need extrapolation
    //T = 800 + 273.15;

    //3. low T, low P
    //T = 50 + 273.15;

    double p_vlh = sw->P_VLH(T);
    double Xl, Xv;
    sw->X_VLH(T, p_vlh, Xl, Xv);
    double Rhol, Hl, dRhodP, dRhodT;
    sw->Rho_phase(T, p_vlh, Xl, Rhol, dRhodP, dRhodT, Liquid);
    sw->H_phase(T, p_vlh, Xl, Hl, Liquid);
}

void debug_VL_props(H2ONaCl::cH2ONaCl* sw)
{
    double T0, P0;

    //1.
    // T0=1.900000e+02 +273.15 , P0=1.254165e+01 * 1E5;

    //2. after Extrapolation, there is still two points difference is very big: T=3.200000e+02 P = 1.127932e+02, Hl=1.461253e+06, see also mmc3 of Driesner2007b. Implement low-T low-p extrapolation, this issue is fixed.
    // T0 = 3.200000e+02+273.14;
    // P0 = 1.127932e+02 * 1E5;

    //3. Rhol big difference: 3.200000e+02 1.127932e+02 6.673641e+02 6.461507e+01 6.027490e+02. This because the salinity in Driesner's result is 0, if I use X=0, will generate the same Rhol.
    // 文件中保存的p值损失精度的原因，导致压力超过了水的boiling压力，所以导致Xl为负，所以导致密度错误，在代码中加入了判断，如果超过量在100Pa以内（这种情况如果是代码内部运算，是自洽的，不会出现问题），则抛出警告；如果超出更多，则直接报错。
    // T0 = 3.200000e+02 + 273.15;
    // P0 = 1.127932e+02 * 1E5;

    //4. Hl big difference: 3.800000e+02,2.345263e+02,8.160492e-04,6.811442e-04,1.349050e-04, ... ,2.121087e+06,2.120716e+06,3.710000e+02 其实质是盐度不同导致的。所以检查这种条件下Xl为何不同,还没找出为什么不同，没法对比。但是与matlab代码对比发现，与matlab代码的结果相同。所以。。。
    // T0 = 3.800000e+02 + 273.15;
    // P0 = 2.345263e+02 * 1E5;

    //5. Rho_v big difference(high-T, low-p at vapor branch): 8.000000e+02,2.464123e+00,2.094581e-04,2.094580e-04,1.000000e-10,1.361844e+02,4.983759e-01,
    T0 = 8.000000e+02 + 273.15;
    P0 = 2.464123e+00 * 1E5;

    double XL = sw->XL_VL(T0,P0);
    double XV = sw->XV_VL(T0,P0);
    double XL_mol = sw->Wt2Mol(XL);
    double XV_mol = sw->Wt2Mol(XV);
    double Rhol,Rhov, Hl,Hv, dRhodP_l, dRhodT_l, dRhodP_v, dRhodT_v;
    sw->H_phase(T0,P0,XL,Hl, Liquid);
    sw->H_phase(T0,P0,XV,Hv, Vapor);
    sw->Rho_phase(T0,P0,XL,Rhol, dRhodP_l, dRhodT_l, Liquid);
    sw->Rho_phase(T0,P0,XV,Rhov, dRhodP_v, dRhodT_v, Vapor);

    // DeformLinearMesh pb_l = sw->PhaseBoundary_VL_DeformLinear(Liquid);
    // DeformLinearMesh pb_v = sw->PhaseBoundary_VL_DeformLinear(Vapor);
    // for (int i = 0; i < pb_l.T.size(); ++i) {
    //     for (int j = 0; j < pb_l.T[i].size(); ++j) {
    //         // double p_crit;
    //         // double p_h2o = sw->get_pWater()->Boiling_p(pb_l.T[i][j]);
    //         // sw->P_Critical(pb_l.T[i][j],p_crit);
    //         // std::cout<<pb_l.T[i][j]<<" "<<pb_l.p[i][j]<<" "<<pb_l.X[i][j]
    //         // <<", P_crit: "<<p_crit
    //         // <<", P_h2o: "<<p_h2o
    //         // <<", T_h2o: "<<sw->get_pWater()->T_critical()
    //         // <<", dp: "<<pb_l.p[i][j]-p_h2o
    //         // <<std::endl;
    //         double Rhol,Hl, Rhov, Hv;
    //         sw->Rho_phase(pb_l.T[i][j], pb_l.p[i][j], pb_l.X[i][j],Rhol, Liquid);
    //         sw->Rho_phase(pb_v.T[i][j], pb_v.p[i][j], pb_v.X[i][j],Rhov, Vapor);
    //         sw->H_phase(pb_l.T[i][j], pb_l.p[i][j], pb_l.X[i][j],Rhol, Liquid);
    //         sw->H_phase(pb_v.T[i][j], pb_v.p[i][j], pb_v.X[i][j],Rhov, Vapor);
    //     }
    //     // std::vector<double> rhol = sw->Rho_phase(pb_l.T[i],pb_l.p[i], pb_l.X[i],Liquid);
    //     // std::vector<double> rhov = sw->Rho_phase(pb_v.T[i],pb_v.p[i], pb_v.X[i],Vapor);
    //
    // }

}

void test_HaliteLiquidus(H2ONaCl::cH2ONaCl* sw)
{
    TriMesh pb_haliteLiquidus = sw->PhaseBoundary_HaliteLiquidus("txt");
}

void test_criticalCurve(H2ONaCl::cH2ONaCl* sw)
{
    STATUS("Properties at the lowest critical point");
    double X, rho, h, dRhodP_l, dRhodT_l;
    sw->Rho_phase(sw->get_pWater()->T_critical(), sw->get_pWater()->p_critical(), 0, rho, dRhodP_l, dRhodT_l, Liquid);
    sw->H_phase(sw->get_pWater()->T_critical(), sw->get_pWater()->p_critical(), 0, h, Liquid);
    cout<<"X: "<<X<<", rho: "<<rho<<", h: "<<h<<endl;
}

void test_constants( H2ONaCl::cH2ONaCl* sw)
{
    cout<<sw->get_pWater()<<endl;
    cout<<"Tmin [K]: "<<sw->Tmin()<<endl;
    cout<<"Tmax [K]: "<<sw->Tmax()<<endl;
    cout<<"pmin [Pa]: "<<sw->pmin()<<endl;
    cout<<"pmax [Pa]: "<<sw->pmax()<<endl;
}