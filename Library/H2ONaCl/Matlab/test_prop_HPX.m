% add xThermal path
addpath('/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/xThermal/Library/API/Matlab');
addpath('/Users/zguo/MyData/Research/3_CodeProject/Hydrothermal-OpenFOAM/H2O_NaCl_EOS_Driesner_2007');
% if use shared compiled mex function, please add libxThermal.dylib,
% libgsl.dylib, libCoolProp.dylib to the dynalic library search path.


% [PROP,props] = SinglePointTest(1, 100E5);
% [Rho, H, Mu, props] = SinglePointTest_water(1000,100E5,'IAPWS95');

% test_constP_H2ONaCl(100E5, 'IAPWS95');
test_constX_H2ONaCl(1e-6, 'IAPS84');
   
% [Rho, H, Mu, props] = test_prop_Water('IAPWS95');


function [Rho, H, Mu, props] = SinglePointTest_water(T,P,name_EOS)
    % T: deg.C, P: Pa, X: kg/kg
    visc_on= 1;
    rho0 = 1;
    [ Rho, dRhodP, H, Mu ] = water_tp_IAPS84(P,T, rho0, 1e-9, 1 ); 

    % xThermal
    [props] = prop_water_TP(T+273.15, P, name_EOS);
end

function [PROP, props] = SinglePointTest(T,P,X)
    % T: deg.C, P: Pa, X: kg/kg
    visc_on= 1;
    [ PROP ] = fluidprop_NaCl_PTX(P,T,X,visc_on);

    % xThermal
    props = prop_TPX(T + 273.15, P, X);
end

function test_constX_H2ONaCl(X0, backendName)
    % calculate properties in pTX space
    visc_on = 1;

    H = linspace(0.1E6,4.5E6,100);
    P = linspace(1,500,250)*1E5;
    [HH,PP] = meshgrid(H,P);
    XX=PP*0 + X0;

    tic
    % calculate
    % [ PROP ] = fluidprop_NaCl_PTX( P,T, Xw, visc_on );
    [ PROP ] = fluidprop_NaCl_PhX(reshape(PP,[],1) , reshape(HH,[],1), reshape(XX,[],1),visc_on);
    toc
    
    % plot
    rows=3;
    cols=5;
    figure(200);clf; 
    
    subplot(rows,cols,1); contourf(HH,PP,reshape(PROP.Reg,size(PP)));colorbar; title('phase region: matlab'); 
    subplot(rows,cols,2); contourf(HH,PP,reshape(PROP.Rho,size(PP)));colorbar; title('Rho: matlab'); 
    subplot(rows,cols,3); contourf(HH,PP,reshape(PROP.h,size(PP)));colorbar; title('H: matlab'); 
    subplot(rows,cols,4); contourf(HH,PP,reshape(PROP.mu_l,size(PP)));colorbar; title('Mu_l: matlab');
    subplot(rows,cols,5); contourf(HH,PP,reshape(PROP.mu_v,size(PP)));colorbar; title('Mu_v: matlab');
    % plot(T, props.Rho,'r','LineWidth',6);
    % plot(T, PROP.Rho,'b-','LineWidth',3);
    
    % use xThermal API
    % [P,X] = P_X_Critical(linspace(0,800,10)+273.15)
    % [P,X] = P_X_Critical(T+273.15);
    % [X_L, X_V] = X_VL(T+273.15, P,backend_name);
    % [Rho_L, Rho_V] = Rho_VL(T+273.15, P, X_L, X_V,backend_name);
    tic
    % props = prop_TPX(T + 273.15, P, Xw);
    try
        props = prop_HPX(HH, PP, XX, backendName);
    catch 
        printf(0, 'error');
    end
    toc

    subplot(rows,cols,6); contourf(HH,PP,props.phase);colorbar; title('phase region: xThermal:PROST'); 
    subplot(rows,cols,7); contourf(HH,PP,props.Rho);colorbar; title('Rho: xThermal:PROST'); 
    subplot(rows,cols,8); contourf(HH,PP,props.H);colorbar; title('H: xThermal:PROST'); 
    subplot(rows,cols,9); contourf(HH,PP,props.Mu_l);colorbar; title('Mu_l: xThermal:PROST'); 
    subplot(rows,cols,10); contourf(HH,PP,props.Mu_v);colorbar; title('Mu_v: xThermal:PROST'); 
    
    % diff
    subplot(rows,cols,11); contourf(HH,PP,(reshape(PROP.Reg,size(PP)) - props.phase));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,12); contourf(HH,PP,(reshape(PROP.Rho,size(PP)) - props.Rho));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,13); contourf(HH,PP,(reshape(PROP.h,size(PP)) - props.H));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,14); contourf(HH,PP,(reshape(PROP.mu_l,size(PP)) - props.Mu_l));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,15); contourf(HH,PP,(reshape(PROP.mu_v,size(PP)) - props.Mu_v));colorbar; title('matlab - xThermal'); 

    hold off;
end

function test_constP_H2ONaCl(P0, backendName)
    % calculate properties in pTX space
    visc_on = 1;

    T = linspace(0.1,900,100);
    % X = 10.^linspace(-10,0,100);
    X = linspace(0.001,1,100);
    [TT,XX] = meshgrid(T,X);
    PP=TT*0 + P0;

    tic
    % calculate
    % [ PROP ] = fluidprop_NaCl_PTX( P,T, Xw, visc_on );
    [ PROP ] = fluidprop_NaCl_PTX(reshape(PP,[],1) , reshape(TT,[],1), reshape(XX,[],1),visc_on);
    toc
    
    % plot
    rows=3;
    cols=5;
    figure(200);clf; 
    
    subplot(rows,cols,1); contourf((XX),TT,reshape(PROP.Reg,size(TT)));colorbar; title('phase region: matlab'); 
    subplot(rows,cols,2); contourf((XX),TT,reshape(PROP.Rho,size(TT)));colorbar; title('Rho: matlab'); 
    subplot(rows,cols,3); contourf((XX),TT,reshape(PROP.h,size(TT)));colorbar; title('H: matlab'); 
    subplot(rows,cols,4); contourf((XX),TT,reshape(PROP.mu_l,size(TT)));colorbar; title('Mu_l: matlab');
    subplot(rows,cols,5); contourf((XX),TT,reshape(PROP.mu_v,size(TT)));colorbar; title('Mu_v: matlab');
    % plot(T, props.Rho,'r','LineWidth',6);
    % plot(T, PROP.Rho,'b-','LineWidth',3);
    
    % use xThermal API
    % [P,X] = P_X_Critical(linspace(0,800,10)+273.15)
    % [P,X] = P_X_Critical(T+273.15);
    % [X_L, X_V] = X_VL(T+273.15, P,backend_name);
    % [Rho_L, Rho_V] = Rho_VL(T+273.15, P, X_L, X_V,backend_name);
    tic
    % props = prop_TPX(T + 273.15, P, Xw);
    props = prop_TPX(TT + 273.15, PP, XX, backendName);
    toc

    subplot(rows,cols,6); contourf((XX),TT,props.phase);colorbar; title('phase region: xThermal:PROST'); 
    subplot(rows,cols,7); contourf((XX),TT,props.Rho);colorbar; title('Rho: xThermal:PROST'); 
    subplot(rows,cols,8); contourf((XX),TT,props.H);colorbar; title('H: xThermal:PROST'); 
    subplot(rows,cols,9); contourf((XX),TT,props.Mu_l);colorbar; title('Mu_l: xThermal:PROST'); 
    subplot(rows,cols,10); contourf((XX),TT,props.Mu_v);colorbar; title('Mu_v: xThermal:PROST'); 
    
    % diff
    subplot(rows,cols,11); contourf((XX),TT,(reshape(PROP.Reg,size(TT)) - props.phase));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,12); contourf((XX),TT,(reshape(PROP.Rho,size(TT)) - props.Rho));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,13); contourf((XX),TT,(reshape(PROP.h,size(TT)) - props.H));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,14); contourf((XX),TT,(reshape(PROP.mu_l,size(TT)) - props.Mu_l));colorbar; title('matlab - xThermal'); 
    subplot(rows,cols,15); contourf((XX),TT,(reshape(PROP.mu_v,size(TT)) - props.Mu_v));colorbar; title('matlab - xThermal'); 

    hold off;
end

function [Rho, H, Mu, props] = test_prop_Water(name_EOS_xThermal)
    T = linspace(0.1,1000,100);
    P = linspace(1,500,100)*1E5;
    [TT,PP] = meshgrid(T,P);
    rho0 = TT*0 + 1;
    [ Rho, dRhodP, H, Mu ] = water_tp_IAPS84(reshape(PP,[],1), reshape(TT,[],1), reshape(rho0, [],1), 1e-9, 1 ); 
    Rho = reshape(Rho, size(TT));
    H = reshape(H, size(TT));
    Mu = reshape(Mu, size(TT));
    
    figure(100);clf;
    rows=3;
    cols=3;
    %plot
    subplot(rows, cols, 1); contourf(TT,PP/1E5,Rho);colorbar; title('Rho: matlab'); 
    subplot(rows, cols, 2); contourf(TT,PP/1E5,H);colorbar; title('H: matlab'); 
    subplot(rows, cols, 3); contourf(TT,PP/1E5,log10(Mu));colorbar; title('Mu: matlab'); 
    fprintf(1,"Mu_min = %.3E, Mu_max = %.3E\n", min(min(Mu)),max(max(Mu)));
    nan_Mu = Mu(isnan(Mu));
    %xThermal
    [props] = prop_water_TP(TT+273.15, PP, name_EOS_xThermal);
    subplot(rows, cols, 4); contourf(TT,PP/1E5,props.Rho);colorbar; title(['Rho: xThermal:',name_EOS_xThermal]); 
    subplot(rows, cols, 5); contourf(TT,PP/1E5,props.H);colorbar; title(['H: xThermal:',name_EOS_xThermal]); 
    subplot(rows, cols, 6); contourf(TT,PP/1E5,log10(props.Mu));colorbar; title(['Mu: xThermal:',name_EOS_xThermal]); 
    fprintf(1,"Mu_min = %.3E, Mu_max = %.3E\n", min(min(props.Mu)),max(max(props.Mu)));
    nan_Mu_PROST = props.Mu(isnan(props.Mu));
    
    % diff
    subplot(rows, cols, 7); contourf(TT,PP/1E5,(Rho - props.Rho)./Rho);colorbar; title('Rho: matlab - xThermal'); 
    subplot(rows, cols, 8); contourf(TT,PP/1E5,(H - props.H)./H);colorbar; title('H: matlab - xThermal'); 
    subplot(rows, cols, 9); contourf(TT,PP/1E5,(Mu - props.Mu)./Mu);colorbar; title('Mu: matlab - xThermal'); 
end