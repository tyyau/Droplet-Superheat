clc; clearvars; close all;

addpath('../include/CorrelationProp');

%% Data
%%
% palette = {'NC16H34','NC7H16'};
% fuel = CorrelationProp(palette);
% load('../Data/hept50hex50_2bar.mat');
% % TODO : Construct based on palette
% W_Comp = [W_NC16H34 W_NC7H16];
% P0 = 2e5;
% SuperheatPlot(fuel, W_Comp, t, r, T, P0, 0)
% %%
% load('Data/hept50hex50_5bar.mat')
% W_Comp = [W_NC16H34 W_NC7H16];
% Comp_Prop = [NC16H34_Prop NC7H16_Prop];
% P0 = 5;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
% %r=0.254
% %%
% load('Data/wat10hex90.mat')
% W_Comp = [W_H2O W_NC16H34];
% Comp_Prop = [H2O_Prop NC16H34_Prop];
% P0 = 6;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)
% %%
% load('Data/eth50dec50_1bar.mat')
% W_Comp = [W_C2H5OH W_NC10H22];
% Comp_Prop = [C2H5OH_Prop NC10H22_Prop];
% P0 = 1;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)
% %%
% load('Data/eth50dec50_10bar.mat')
% W_Comp = [W_C2H5OH W_NC10H22];
% Comp_Prop = [C2H5OH_Prop NC10H22_Prop];
% P0 = 10;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)
% %%
% load('Data/meth50dod50.mat')
% W_Comp = [W_NC12H26 W_CH3OH];
% Comp_Prop = [NC12H26_Prop CH3OH_Prop];
% P0 = 40;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
% %%
% load('Data/meth50hex50.mat')
% W_Comp = [W_NC16H34 W_CH3OH];
% Comp_Prop = [NC16H34_Prop CH3OH_Prop];
% P0 = 40;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
% %%
% load('Data/meth40hex60.mat')
% W_Comp = [W_NC16H34 W_CH3OH];
% Comp_Prop = [NC16H34_Prop CH3OH_Prop];
% P0 = 40;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
% %%
% load('Data/hept50hex50_4bar.mat')
% W_Comp = [W_NC16H34 W_NC7H16];
% Comp_Prop = [NC16H34_Prop NC7H16_Prop];
% P0 = 4;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
% %r=0.21
% %%
% load('Data/hept50hex50_3bar.mat')
% W_Comp = [W_NC16H34 W_NC7H16];
% Comp_Prop = [NC16H34_Prop NC7H16_Prop];
% P0 = 3;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
% %%
% load('Data/eth50dod50_1bar.mat')
% W_Comp = [W_C2H5OH W_NC12H26];
% Comp_Prop = [C2H5OH_Prop NC12H26_Prop];
% P0 = 1;
% SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)

%% Kinetic Model
% load('../Data/hept50hex50_5bar.mat')
% palette = {'NC16H34','NC7H16'};
% fuel = CorrelationProp(palette);
% W_Comp = [W_NC16H34 W_NC7H16];
% P0 = 5e5;
% KineticModel(fuel, W_Comp, t, r, T, P0, 0)

%%
load('../Data/eth50dec50_1bar.mat')
palette = {'C2H5OH','NC10H22'};
fuel = CorrelationProp(palette);
W_Comp = [W_C2H5OH W_NC10H22];
P0 = 1e5;
KineticModel(fuel, W_Comp, t, r, T, P0, 0)
%%
load('Data/hept50hex50_2bar.mat')
W_Comp = [W_NC16H34 W_NC7H16];
Comp_Prop = [NC16H34_Prop NC7H16_Prop];
Antoine = [Ant_NC16H34 Ant_NC7H16];
P0 = 2;
KineticModel(fuel, Antoine, W_Comp, t, r, T, P0, 0)
%%
load('Data/water_9K.mat')
W_Comp = [W_H2O];
Comp_Prop = [H2O_Prop];
Antoine = [Ant_H2O];
P0 = .3866;
[~,~,tf1,rf1]=KineticModel(Comp_Prop, Antoine, W_Comp, t, r, T, P0, 0);

load('Data/water_11K.mat')
W_Comp = [W_H2O];
Comp_Prop = [H2O_Prop];
Antoine = [Ant_H2O];
P0 = .1266;
[~,~,tf2,rf2]=KineticModel(Comp_Prop, Antoine, W_Comp, t, r, T, P0, 0);

load('Data/water_droplet_data.mat')

plot(tf1,rf1*10,tf2,rf2*10,t_circle, r_circle, '*', t_square, r_square, '*',t_circlez, r_circlez, '--',t_squarez, r_squarez, '--');
ylabel('radius[cm]');
xlabel('time[s]');
legend('P = 38kPa', 'P = 13kPa', 'P = 38kPa', 'P = 13kPa');
title('Water bubble growth verification');