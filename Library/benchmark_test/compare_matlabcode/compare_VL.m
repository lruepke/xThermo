% calculate phase boundary of VL and related density & enthalpy using function calc_boundary_prop.m of the matlab code : https://github.com/lruepke/H2O_NaCl_EOS_Driesner_2007

PhaseBoundary_VL();

function PhaseBoundary_VL()
    Tcrit_h2o = 373.976;
    T1 = linspace(0.1,Tcrit_h2o,200);
    T2 = linspace(Tcrit_h2o, 1000, 200);
    T = cat(2, T1, T2);
    for i=1:length(T)
        P_crit(i) = p_crit_H2ONaCl(T(i));
        P_vlh(i) = calc_P_vlh(T(i));
        if(T(i)<Tcrit_h2o)
            P_h2o = fluidprop_crit_T( T(i), 1e-8) * 1E6; %fluidprop_crit_T return pressure in MPa
            P_crit(i) = P_h2o;
        end
        P_vlh(P_vlh<=0)=0.1E5;
    end
    
    fpout_T=fopen('TT_VL.txt','w');
    fpout_p=fopen('PP_VL.txt','w');
    fpout_Xl=fopen('Xl_VL.txt','w');
    fpout_Xv=fopen('Xv_VL.txt','w');
    for iT = 1:length(T)
        P = linspace(P_vlh(iT),P_crit(iT),200);
        for ip = 1:length(P)
            [  Xw_l, Xw_v ] = calc_boundary_prop( P(ip), T(iT), 1 );
            fprintf(fpout_T,'%.8E ', T(iT)+273.15);
            fprintf(fpout_p,'%.8E ', P(ip));
            fprintf(fpout_Xl,'%.8E ', Xw_l);
            fprintf(fpout_Xv,'%.8E ', Xw_v);
        end
        fprintf(fpout_T,'\n');
        fprintf(fpout_p,'\n');
        fprintf(fpout_Xl,'\n');
        fprintf(fpout_Xv,'\n');
    end
    fclose(fpout_T);
    fclose(fpout_p);
    fclose(fpout_Xl);
    fclose(fpout_Xv);
end