clc; clear all; close all;

% ADD PATHS
addpath('../include/CorrelationProp');

% PALETTE AND INITIAL LIST MUST BE INPUT FOR CORRELATION MODE
% palette = funcName;
% initList = initName;
palette = {'IC8H18','C7H8','C7H16'};
initList = [0.56,0.28,0.16];
p = 101325; % Pa
T = 300; % K

fuel = CorrelationProp(palette,initList);

% PRINT ALL PROPERTIES
disp(palette);
fprintf('MW = %s\n',sprintf('%.4g ',fuel.MW));
fprintf('Tc = %s\n',sprintf('%.4g ',fuel.TcVec));
fprintf('Pc = %s\n',sprintf('%.4g ',fuel.PcVec));
fprintf('Acentric = %s\n',sprintf('%.4g ',fuel.omegaVec));
fprintf('Latent Heat = %s\n',sprintf('%.4g ',fuel.L(T)));
fprintf('Specific Heat = %s\n',sprintf('%.4g ',fuel.c_l(T)));
fprintf('Diffusivity = %s\n',sprintf('%.4g ',fuel.D(p,T)));
fprintf('Saturated Vapor Pressure = %s\n',sprintf('%.4g ',fuel.Psat(T)));
fprintf('Thermal Conductivity = %s\n',sprintf('%.4g ',fuel.lambdaL(T)));
fprintf('Surface Tension = %s\n',sprintf('%.4g ',fuel.sigma(T)));

