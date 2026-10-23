clear; clc; clf;
% ==============

file_LUT_TPX = 'lut_constX_TP_9.bin';
% file_LUT_TPX = 'TPX/lut_TPX_9.bin';
file_LUT_HPX = 'lut_constX_HP_10.bin';
% file_LUT_HPX = 'HPX/lut_HPX_9.bin';
file_LUT_HP_water = 'lut_constX_HP_11.bin';

[H,P,X, props] = test_lookup_HPX(file_LUT_HPX);
% [T,P,X, props] = test_lookup_TPX(file_LUT_TPX);
% water
% [H,P,props] = test_lookup_HP(file_LUT_HP_water);

function [H,P,props] = test_lookup_HP(file_LUT)
    fluidName = 'H2O';
    H0 = 1E6; % J/kg
    P0 = 100E5; %Pa
%     X0 = 0.1; %kg/kg
    show = true;
    % 1. single point lookup 
    H=H0;
    P=P0;
    props = lookup(file_LUT, fluidName, H0,P0);

    % 2. vector lookup
    H = linspace(0.1E6, 4E6, 100);
    P = H*0 + P0;
    props = lookup(file_LUT, fluidName, H,P);
    if(show)
        figure(2);
        plot(H, props.Rho);
    end

    % 3. matrix lookup
    H = linspace(0.1E6, 4E6, 100);
    P = linspace(10E5, 500E5, 100);
    [HH,PP] = meshgrid(H,P);
    props = lookup(file_LUT, fluidName, HH,PP);
    prop_names = fieldnames(props);
    if(show)
        figure(3);
        ncol = 2;
        nrow = ceil(length(prop_names)/ncol);
        for i = 1:length(prop_names)
            subplot(nrow, ncol, i);
            hold on;
            pcolor(HH/1E6,PP/1E5,props.(cell2mat(prop_names(i))));
            shading interp;
            if(~strcmp(cell2mat(prop_names(i)) , 'Phase'))
                contour(HH/1E6,PP/1E5,props.(cell2mat(prop_names(i))), 'LineColor', 'k');
            end
            hold off;
            title(prop_names(i));
            xlabel('Bulk specific enthalpy (MJ/kg)');
            ylabel('Pressure (bar)');
        end
    end
end

function [H,P,X, props] = test_lookup_HPX(file_LUT)
    fluidName='H2O-NaCl';
    H0 = 1E6; % J/kg
    P0 = 100E5; %Pa
    X0 = 0.1; %kg/kg
    show = true;
    % 1. single point lookup 
    H=H0;
    P=P0;
    X=X0;
    props = lookup(file_LUT,fluidName, H0,P0,X0);

    % 2. vector lookup
    H = linspace(0.1E6, 4E6, 100);
    P = H*0 + P0;
    X = H*0 + X0;
    props = lookup(file_LUT,fluidName, H,P,X);
    if(show)
        figure(2);
        plot(H, props.Rho);
    end

    % 3. matrix lookup
    H = linspace(0.1E6, 4E6, 100);
    P = linspace(10E5, 500E5, 100);
    [HH,PP] = meshgrid(H,P);
    XX = HH*0 + X0;
    props = lookup(file_LUT, fluidName, HH,PP,XX);
    prop_names = fieldnames(props);
    if(show)
        figure(3);
        ncol = 2;
        nrow = ceil(length(prop_names)/ncol);
        for i = 1:length(prop_names)
            subplot(nrow, ncol, i);
            hold on;
            pcolor(HH/1E6,PP/1E5,props.(cell2mat(prop_names(i))));
            shading interp;
            if(~strcmp(cell2mat(prop_names(i)) , 'Phase'))
                contour(HH/1E6,PP/1E5,props.(cell2mat(prop_names(i))), 'LineColor', 'k');
            end
            hold off;
            title(prop_names(i));
            xlabel('Bulk specific enthalpy (MJ/kg)');
            ylabel('Pressure (bar)');
        end
    end
end

function [T,P,X, props] = test_lookup_TPX(file_LUT)
    fluidName = 'H2O-NaCl';
    T0 = 100 + 273.15;
    P0 = 100E5;
    X0 = 0.1;
    show = true;
    % 1. single point lookup 
    props0 = lookup(file_LUT, fluidName, T0,P0,X0);

    % 2. vector lookup
    T = linspace(1, 500, 100) + 273.15;
    P = T*0 + P0;
    X = T*0 + X0;
    props = lookup(file_LUT, fluidName, T,P,X);
    if(show)
        figure(2);
        plot(T, props.Rho);
    end

    % 3. matrix lookup
    T = linspace(1, 500, 100) + 273.15;
    P = linspace(10E5, 500E5, 100);
    [TT,PP] = meshgrid(T,P);
    XX = TT*0 + X0;
    props = lookup(file_LUT, fluidName, TT,PP,XX);
    prop_names = fieldnames(props);
    if(show)
        figure(3);
        ncol = 2;
        nrow = ceil(length(prop_names)/ncol);
        for i = 1:length(prop_names)
            subplot(nrow, ncol, i);
            hold on;
            pcolor(TT-273.15,PP/1E5,props.(cell2mat(prop_names(i))));
            shading interp;
            if(~strcmp(cell2mat(prop_names(i)) , 'Phase'))
                contour(TT-273.15,PP/1E5,props.(cell2mat(prop_names(i))), 'LineColor', 'k');
            end
            hold off;
            title(prop_names(i));
            xlabel('Temperature(deg.C)');
            ylabel('Pressure (bar)');
        end
    end
end