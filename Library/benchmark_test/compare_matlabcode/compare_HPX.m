
constX(0.2);

function constX(x)
    tol = 1E-10;
    Hmin = 0.1E6;
    Hmax = 4.5E6;
    Pmin = 1E5;
    Pmax = 600E5;

    H = linspace(Hmin, Hmax, 100);
    P = linspace(Pmin, Pmax ,100);
    fpout_H=fopen('HH.txt','w');
    fpout_p=fopen('PP.txt','w');
    fpout_Rho=fopen('RHO.txt','w');
    fpout_T=fopen('TT.txt','w');
    for i = 1:length(H)
        h = H(i);
        for j = 1:length(P)
            p = P(j);
            [ PROP ] = fluidprop_NaCl_PhX( p, h, x,tol );
            fprintf(fpout_H,'%.8E ', H);
            fprintf(fpout_p,'%.8E ', p);
            fprintf(fpout_Rho,'%.8E ', PROP.rho);
            fprintf(fpout_T,'%.8E ', PROP.T);
        end
        fprintf(fpout_H,'\n');
        fprintf(fpout_p,'\n');
        fprintf(fpout_Rho,'\n');
        fprintf(fpout_T,'\n');
    end
    fclose(fpout_H);
    fclose(fpout_p);
    fclose(fpout_Rho);
    fclose(fpout_T);
end