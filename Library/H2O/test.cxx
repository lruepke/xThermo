/**
 * @file test.cpp
 * @author Zhikui Guo (zguo@geomar.de)
 * @brief Test implementation of H2O
 * @version 0.1
 * @date 2022-03-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */
// #include "thermo.h"
#include "IAPS84.h"
#include "IAPWS95_CoolProp.h"
#include "IAPWS95.h"
#include "IAPWS-IF97.h"
#include "LookUpTableForestI.H"

using namespace std;
using namespace xThermal;

void test_constants( xThermal::cxThermal* thermo);
void test_PhaseDiagram(xThermal::cxThermal* thermo, double T, double p);
void test_Rho(xThermal::cxThermal& thermo, const double& T, const double& P, double& Rho);
void test_Performance();
void test_water_profile(xThermal::cxThermal* thermo, int nT=100, int nP = 100);
int test_LUTgen(int argc, char** argv, xThermal::cxThermal& thermo);
int lutinfo(int argc, char** argv);
int test_LUT_lookup(int argc, char** argv, xThermal::cxThermal& thermo);
void test_viscosity(xThermal::cxThermal* thermo);
void test_PH(xThermal::cxThermal* thermo);
void test_compressibility_twophase_in_porousmedium(xThermal::cxThermal* thermo, double S=0.5, double porosity=0.1, double rho_r=2000, double cp_r=1000);
int main(int argc, char** argv)
{
    PROST::cIAPS84 iaps84;
    IAPWS95::cIAPWS95 iapws95;
#ifdef USE_COOLPROP
    COOLPROP::cIAPWS95_CoolProp iapws95_coolprop;
#endif

    // iapws95.Verification_Mu();
    // test_viscosity(&iapws95);
    // test_constants(&iaps84);
    // test_PhaseDiagram(&iaps84, 100+273.15, 200E5);
    // test_PhaseDiagram(&iaps84, 600+273.15, 200E5);
    // test_PhaseDiagram(&iaps84, iaps84.T_critical(), iaps84.p_critical());
    // iapws95.loadLUT("lut_constX_TP_12.bin");
    // test_LUTgen(argc, argv,iaps84);
    // lutinfo(argc,argv);
    //    test_LUT_lookup(argc,argv, iaps84);
    //    exit(0);


    // test_constants(&iapws95);
    //  test_PhaseDiagram(&iapws95_coolprop, 100+273.15, 200E5);
    //  test_PhaseDiagram(&iapws95_coolprop, 600+273.15, 200E5);
    //  test_PhaseDiagram(&iapws95_coolprop, iapws95.T_critical(), iapws95.p_critical());
    //    test_LUTgen(argc, argv,iapws95_coolprop);

    //    ThermodynamicProperties props;
    ////    iapws95_coolprop.UpdateState_TPX(props,554.936,iapws95_coolprop.Boiling_p(554.936));
    //    cout<<iapws95_coolprop.phase_name(iapws95_coolprop.findPhaseRegion_TPX(554.936,101E5))<<endl;
    //    exit(0);

    //    ThermodynamicProperties props;
    //    double p = CoolProp::saturation_ancillary("Water","P", 1, "T", T0);
    //    cout<<CoolProp::PhaseSI("T", T0,"P", p, "Water")<<endl;
    //    iapws95_coolprop.UpdateState_TPX(props, T0, p_sat);
    //    cout<<"Rho 95: "<<props.Rho<<" p: "<<p<<" boil p: "<<iapws95_coolprop.Boiling_p(T0)<<endl;


    // test_constants(&iapws95);
    // test_PhaseDiagram(&iapws95, 100+273.15, 200E5);

    // IAPWS_IF97::cIAPWS_IF97 iapws_if97;
    // test_constants(&iapws_if97);
    // test_PhaseDiagram(&iapws_if97, 100+273.15, 200E5);
    //     double Rho, Rho_min, Rho_max;
    //     iapws95.loadLUT("lut_constX_TP_12.bin");
    //     test_water_profile(&iaps84);
    //     test_water_profile(&iapws95);
    // //    test_water_profile(iapws95_coolprop);
    // test_PH(&iaps84);
    // test_PH(&iapws95);
    // test_PH(&iapws95_coolprop);
    // ThermodynamicProperties props;
    // double h0=518008.007503, p0=10002936.627387;
    // iapws95.UpdateState_HPX(props,h0, p0);
    // std::cout<<props<<endl;
    // iaps84.UpdateState_HPX(props,h0, p0);
    // std::cout<<props<<endl;
    // props = iapws95.UpdateState_TPX(532.163, 3.99239e+06);
    // std::cout<<props<<endl;

    // test_compressibility_twophase_in_porousmedium(&iapws95_coolprop);

    return 0;
}

void test_compressibility_twophase_in_porousmedium(xThermal::cxThermal* thermo, double S, double porosity, double rho_r, double cp_r)
{
    // the default parameters are taken from Grant et al.(1979, doi: 10.1029/wr015i003p00684 ), below the Eq.(7)
    double T0 = 250 + 273.15;
    ThermodynamicProperties props_boil;
    thermo->Boiling_p(T0, props_boil);
    double L = props_boil.H_v - props_boil.H_l; // Latent heat is needed for phase changes from liquid to vapor
    double dT = 0.1;
    double dpdT = (thermo->Boiling_p(T0+dT/2) - thermo->Boiling_p(T0-dT/2))/dT;
    double compressibility = ((1-porosity)*rho_r*cp_r + porosity*S*props_boil.Rho_l*props_boil.Cp_l) * (props_boil.Rho_l - props_boil.Rho_v) / (porosity*L*dpdT*props_boil.Rho_l*props_boil.Rho_v); // Grant(1979) formula, Eq.(7)
    double compressibility_withVapor = ((1-porosity)*rho_r*cp_r + porosity*S*props_boil.Rho_l*props_boil.Cp_l + porosity*(1-S)*props_boil.Rho_v*props_boil.Cp_v ) * (props_boil.Rho_l - props_boil.Rho_v) / (porosity*L*dpdT*props_boil.Rho_l*props_boil.Rho_v); // Grant(1979) formula, Eq.(7)
    double compressibility_Geiger = ((1-porosity)*rho_r*cp_r + porosity*S*props_boil.Rho_l*props_boil.Cp_l) / porosity * pow((props_boil.Rho_l - props_boil.Rho_v)/(L*props_boil.Rho_v*props_boil.Rho_l), 2.0) * T0; // Geiger(2006,doi: 10.1007/s11242-005-0108-z ) Eq.5
    printf("Liquid phase compressibility: %.2E per. bar\n", props_boil.IsothermalCompressibility_l*1E5);
    printf("Vapor phase compressibility: %.2E per. bar\n", props_boil.IsothermalCompressibility_v*1E5);
    printf("Two phase compressibility(Grant): %.3f per. bar\n", compressibility*1E5);
    printf("Two phase compressibility(Grant, with vapor): %.3f per. bar\n", compressibility_withVapor*1E5);
    printf("Two phase compressibility(Geiger): %.3f per. bar\n", compressibility_Geiger*1E5);
    printf("Liquid Cp: %.2E J/kg/K\n", props_boil.Cp_l);
    printf("Vapor Cp: %.2E J/kg/K\n", props_boil.Cp_v);
    printf("dpdT along boiling curve (numerical diff): %.2E Pa/K\n", dpdT);

    cout<<props_boil<<endl;

}

void test_PH(xThermal::cxThermal* thermo)
{
    double P = 100E5, H = 2E6;
    ThermodynamicProperties props = thermo->UpdateState_HPX(H,P);
    cout<<props<<endl;
}
void test_viscosity(xThermal::cxThermal* water)
{
    ThermodynamicProperties props;
    double T,P;
    // T= 990+273.15, P = 300E5;
    // T = 0 + 273.15; P = 100E5;
    // T = 373.470949 + 273.15; P = 219.236019E5; //bug point close Critical point
    T = water->T_critical(), P = water->p_critical();
    water->UpdateState_TPX(props, T,P);
    double mu = props.Mu;
    double Rho = props.Rho;
    double phase = props.phase;
    cout<<props<<endl;
}

int test_LUT_lookup(int argc, char** argv, xThermal::cxThermal& thermo)
{
    using namespace std;
    thermo.loadLUT("lut_constX_TP_10.bin");
//    //simple lookup
//    double* props = new double[((LOOKUPTABLE_FOREST::LookUpTableForest_2D*)thermo.get_pLUT())->m_map_props.size()];
//    double xyz_min_target[2];
//    double T0=573.15, p0=300E5;
//    // 1. lookup
//    thermo.lookup(props, xyz_min_target,T0, p0);
//    cout<<"lut T "<<T0<<", p "<<p0<<", rho: "<<props[0]<<endl;
//    // 2. calculate
//    ThermodynamicProperties prop;
//    thermo.UpdateState_TPX(prop, T0, p0);
//    cout<<"calc T "<<T0<<", p "<<p0<<", rho: "<<prop.Rho<<endl;
//
//    exit(0);

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

    return 0;
}

int lutinfo(int argc, char** argv)
{
    if(argc!=2)ERROR("Usage: lutinfo myLUT.bin");
    // H2ONaCl::cH2ONaCl sw;
    // sw.loadLUT(argv[1]);
    int m_dim_lut = LOOKUPTABLE_FOREST::get_dim_from_binary(argv[1]);

    switch (m_dim_lut)
    {
        case 2:
        {
            LOOKUPTABLE_FOREST::LookUpTableForest_2D lut_2d(argv[1], NULL);
        }
            break;
        case 3:
        {
            LOOKUPTABLE_FOREST::LookUpTableForest_3D lut_3d(argv[1], NULL);
        }
            break;
        default:
        ERROR("The dim in the binary file is neither 2 nor 3, it is not a valid LUT file: "+string(argv[1]));
            break;
    }
    return 0;
}

int test_LUTgen(int argc, char** argv, xThermal::cxThermal& eos)
{
    if(argc!=9)
    {
        STATUS("Usage: "+string(argv[0])+" [T|H] [xmin] [xmax] [ymin] [ymax] [min_level] [max_level] [nThread]");
        STATUS("Unit: T[deg.C], P[bar], X[wt. NaCl, 0-1], H[MJ/kg]");
        STATUS_color("Example 1: T 1 1000 5      600  4   7   8 (var TP)", COLOR_BLUE);
        STATUS_color("Example 5: H 0.1 3.9 5     600  4   7   8 (constX, var HP)", COLOR_PURPLE);
        return 0;
    }
    string TorH = string(argv[1]);
    double xmin = atof(argv[2]);
    double xmax = atof(argv[3]);
    double ymin = atof(argv[4]);
    double ymax = atof(argv[5]);
    int min_level = atoi(argv[6]);
    int max_level = atoi(argv[7]);
    int n_threads = atoi(argv[8]);
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
    double constZ = 0;
    if (TorH == "T")
    {
        eos.createLUT_2D(xmin + 273.15, xmax + 273.15, ymin * 1E5, ymax * 1E5, constZ, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_T, min_level, max_level, Update_prop_CommonlyUsedTPX);
        eos.save_lut_to_binary("lut_constX_TP_"+std::to_string(max_level)+".bin");
        eos.save_lut_to_vtk("lut_constX_TP_"+std::to_string(max_level)+".vtu");
    }else if (TorH == "H")
    {
        eos.createLUT_2D(xmin*1E6, xmax*1E6, ymin*1E5, ymax*1E5, constZ, LOOKUPTABLE_FOREST::CONST_X_VAR_TorHP, LOOKUPTABLE_FOREST::EOS_ENERGY_H, min_level, max_level, Update_prop_CommonlyUsedTPX);
        eos.save_lut_to_binary("lut_constX_HP_"+std::to_string(max_level)+".bin");
        eos.save_lut_to_vtk("lut_constX_HP_"+std::to_string(max_level)+".vtu");
    }else
    {
        ERROR("The first arg must be T or H.");
    }
    return 0;
}

void test_water_profile(xThermal::cxThermal* water, int nT, int nP)
{
    double Tmin = 2 + 273.15, Tmax = 1000 + 273.15;
    double dT = (Tmax - Tmin)/(nT-1);
    double Pmin = 1E5, Pmax = 2000E5;
    double dP = (Pmax - Pmin)/(nT-1);
    time_t start, end;
    ThermodynamicProperties props;
    start = clock();
    for (int i = 0; i < nT; ++i) {
        for (int j = 0; j < nP; ++j) {
            water->UpdateState_TPX(props, Tmin + dT*i, Pmin + dP*j);
        }
    }
    end = clock();
    STATUS_time("Computing time of "+water->name(),end-start);
}
void test_Performance()
{
    int nT = 100000;
    double Tmin = 2 + 273.15, Tmax = 1000 + 273.15;
    double dT = (Tmax - Tmin)/(nT-1);
    double P = 5.057798e+02 * 1E5;
    P = 5E5;
    ThermodynamicProperties props;
    // PROST
    PROST::cIAPS84 iaps84;
    time_t start = clock();
    for (int i = 0; i < nT; ++i) {
        iaps84.UpdateState_TPX(props, Tmin + dT*i,P);
    }
    time_t end = clock();
    STATUS_time("Computing time of "+iaps84.name(),end-start);

    // CoolProp
#ifdef USE_COOLPROP
    COOLPROP::cIAPWS95_CoolProp coolProp;
    start = clock();
    for (int i = 0; i < nT; ++i) {
        coolProp.UpdateState_TPX(props, Tmin + dT*i,P);
        
    }
    end = clock();
    STATUS_time("Computing time of "+coolProp.name(),end-start);
#endif
    // // IAPWS95
    // IAPWS95::cIAPWS95 iapws95;
    // start = clock();
    // for (int i = 0; i < nT; ++i) {
    //     iapws95.UpdateState_TPX(Tmin + dT*i,P, props);
    //     // iapws95.rhomass();
    // }
    // end = clock();
    // STATUS_time("Computing time of "+iapws95.name(),end-start);
}


void test_PhaseDiagram(xThermal::cxThermal* thermo, double T, double p)
{
    ThermodynamicProperties props;
    STATUS("Test phase diagram of "+thermo->name());
    cout<<"T = "<<T<<" [K], p = "<<p<<" [Pa]"<<endl;
    thermo->UpdateState_TPX(props, T, p);
    cout<<"phase: "<<props.phase<<", name: "<<thermo->phase_name(props.phase)<<endl;
    cout<<"Boiling p: "<<thermo->Boiling_p(T)<<endl;
}

void test_constants( xThermal::cxThermal* thermo)
{
    STATUS("Test " + thermo->name());
    cout<<"Tmin [K]: "<<thermo->Tmin()<<endl;
    cout<<"Tmax [K]: "<<thermo->Tmax()<<endl;
    cout<<"pmin [Pa]: "<<thermo->pmin()<<endl;
    cout<<"pmax [Pa]: "<<thermo->pmax()<<endl;
    cout<<"T_triple [K]: "<<thermo->Ttriple()<<endl;
    cout<<"T_crit [K]: "<<thermo->T_critical()<<endl;
    cout<<"p_crit [Pa]: "<<thermo->p_critical()<<endl;
    cout<<"rho_crit [kg/m^3]: "<<thermo->rhomass_critical()<<endl;
    cout<<"rho_crit [mol/m^3]: "<<thermo->rhomolar_critical()<<endl;
    cout<<"Molar mass [kg/mol]: "<<thermo->molar_mass()<<endl;
}