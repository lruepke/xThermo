
// calculate water phase diagram data

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

#include "IAPWS-97.H"
// #define IAPWS97_TCRIT 647.096 /* K */
// #define IAPWS97_PCRIT 22.064e6 /* Pa */
// #define IAPWS97_RHOCRIT 322. /* kg/m³ */

// 1. calculate rho, H, saturation in P-T space
// p: Pa, T: C
void calTable_PT(double plim[2], double Tlim[2], double dp=1e5, double dT=1);
void calTable_PH(double plim[2], double Hlim[2], double dp=1e5, double dH=0.1e5);
void singletest()
{
    double T, p, h, rho;
	// SteamState S = freesteam_set_pT(100e5,200+273.15);
    SteamState S = freesteam_set_ph(350e5, 4.5e6);
	rho=freesteam_rho(S);
    T=freesteam_T(S);
    h=freesteam_h(S);
    p=freesteam_p(S);
	printf("p: %lf,  T: %lf   h: %lf   rho: %lf\n",p, T, h, rho);

}
int main()
{
    singletest();

    // // 1. p-T space
    // double plim_Pa[2]={1e5,400e5}, Tlim_C[2]={1e-5,1000};
    // calTable_PT(plim_Pa,Tlim_C);

    // // 2. p-H space
    // double Hlim[2]={0.05e6, 4.5e6};
    // calTable_PH(plim_Pa,Hlim);

    return 0;
}

void calTable_PH(double plim[2], double Hlim[2], double dp, double dH)
{
    string datapath="../../../data/PhaseDiagram/water/pH/";
    FILE* fp_T=fopen((datapath+"TT.dat").c_str(),"w");
	FILE* fp_P=fopen((datapath+"PP.dat").c_str(),"w");
	FILE* fp_rho=fopen((datapath+"RHO.dat").c_str(),"w");
	FILE* fp_X=fopen((datapath+"X.dat").c_str(),"w");
	FILE* fp_cp=fopen((datapath+"CP.dat").c_str(),"w");
	FILE* fp_mu=fopen((datapath+"MU.dat").c_str(),"w");
	// FILE* fp_alpha=fopen("ALPHA.dat","w");
	// FILE* fp_beta=fopen("BETA.dat","w");
	FILE* fp_region=fopen((datapath+"REGION.dat").c_str(),"w");
	FILE* fp_h=fopen((datapath+"H.dat").c_str(),"w");
    double T=0;
    double rho=0, h, region,X,cp,mu;
    for (double p = plim[0]; p < plim[1]; p=p+dp)
    {
        for (double h = Hlim[0]; h < Hlim[1]; h=h+dH)
        {
            SteamState S = freesteam_set_ph(p,h);
            rho=freesteam_rho(S);
            T=freesteam_T(S);
			region=freesteam_region(S);
			X=freesteam_x(S);
			cp=freesteam_cp(S);
			mu=freesteam_mu(S);
            // write to file
			fprintf(fp_T,"%lf ",T - 273.15);
			fprintf(fp_P,"%lf ",p);
			fprintf(fp_rho,"%lf ",rho);
			fprintf(fp_X,"%lf ",X);
			fprintf(fp_cp,"%lf ",cp);
			fprintf(fp_mu,"%lf ",mu);
			// fprintf(fp_alpha,"%.6E ",alpha);
			// fprintf(fp_beta,"%.6E ",beta);
			fprintf(fp_region,"%.2f ",region);
			fprintf(fp_h,"%lf ",h);
        }
        fprintf(fp_T,"\n");
		fprintf(fp_P,"\n");
		fprintf(fp_rho,"\n");
		fprintf(fp_X,"\n");
		fprintf(fp_cp,"\n");
		fprintf(fp_mu,"\n");
		// fprintf(fp_alpha,"\n");
		// fprintf(fp_beta,"\n");
		fprintf(fp_region,"\n");
		fprintf(fp_h,"\n");
    }
    fclose(fp_T);
	fclose(fp_P);
	fclose(fp_rho);
	fclose(fp_X);
	fclose(fp_cp);
	fclose(fp_mu);
	// fclose(fp_alpha);
	// fclose(fp_beta);
	fclose(fp_region);
	fclose(fp_h);
}


void calTable_PT(double plim[2], double Tlim[2], double dp, double dT)
{
    string datapath=".";
    FILE* fp_T=fopen((datapath+"TT.dat").c_str(),"w");
	FILE* fp_P=fopen((datapath+"PP.dat").c_str(),"w");
	FILE* fp_rho=fopen((datapath+"RHO.dat").c_str(),"w");
	// FILE* fp_cp=fopen("CP.dat","w");
	// FILE* fp_mu=fopen("MU.dat","w");
	// FILE* fp_alpha=fopen("ALPHA.dat","w");
	// FILE* fp_beta=fopen("BETA.dat","w");
	FILE* fp_region=fopen((datapath+"REGION.dat").c_str(),"w");
	FILE* fp_h=fopen((datapath+"H.dat").c_str(),"w");
    double T=0;
    double rho=0, h,region;
    for (double p = plim[0]; p < plim[1]; p=p+dp)
    {
        for (double T_C = Tlim[0]; T_C < Tlim[1]; T_C=T_C+dT)
        {
            T = T_C + 273.15;
            SteamState S = freesteam_set_pT(p,T);
            rho=freesteam_rho(S);
            h=freesteam_h(S);
			region=freesteam_region(S);

            // write to file
			fprintf(fp_T,"%lf ",T_C);
			fprintf(fp_P,"%lf ",p);
			fprintf(fp_rho,"%lf ",rho);
			// fprintf(fp_cp,"%lf ",cp);
			// fprintf(fp_mu,"%lf ",mu);
			// fprintf(fp_alpha,"%.6E ",alpha);
			// fprintf(fp_beta,"%.6E ",beta);
			fprintf(fp_region,"%.2f ",region);
			fprintf(fp_h,"%lf ",h);
        }
        fprintf(fp_T,"\n");
		fprintf(fp_P,"\n");
		fprintf(fp_rho,"\n");
		// fprintf(fp_cp,"\n");
		// fprintf(fp_mu,"\n");
		// fprintf(fp_alpha,"\n");
		// fprintf(fp_beta,"\n");
		fprintf(fp_region,"\n");
		fprintf(fp_h,"\n");
    }
    fclose(fp_T);
	fclose(fp_P);
	fclose(fp_rho);
	// fclose(fp_cp);
	// fclose(fp_mu);
	// fclose(fp_alpha);
	// fclose(fp_beta);
	fclose(fp_region);
	fclose(fp_h);
}
