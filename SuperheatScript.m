%% Properties
% IC8H18 = isooctane
% TMBENZ = 1,3,5-trimethylbenzene
% NPBENZ = n-propylbenzene
% NC12H26 = n-dodecane
% NC7H16 = n-heptane
% NC10H22 = n-decane
% C2H5OH = ethanol
% NC16H26 = n-hexadecane
% CH3OH = methanol
% H2O = water

% Perry's properties
% Comp_prop = [M_Comp Tc_Comp Pc_Comp Vc_Comp]';
M_IC8H18 = 114.231;
Tc_IC8H18 = 543.96; % K
Pc_IC8H18 = 25.6; % bar
Vc_IC8H18 = .465; % l/mol
om_IC8H18 = .264;
Tb_IC8H18 = 372.5;
Ant_IC8H18 = [6.99021 1358.75 231.405]'; % Knovel
IC8H18_Prop = [M_IC8H18 Tc_IC8H18 Pc_IC8H18 Vc_IC8H18 om_IC8H18 Tb_IC8H18]';
M_TMBENZ = 120.194;
Tc_TMBENZ = 637.36;
Pc_TMBENZ = 31.1;
Vc_TMBENZ = .433;
om_TMBENZ = .397;
Tb_TMBENZ = 437.8;
Ant_TMBENZ = [7.26105 1695.83 222.415];
TMBENZ_Prop = [M_TMBENZ Tc_TMBENZ Pc_TMBENZ Vc_TMBENZ om_TMBENZ Tb_TMBENZ]';
M_NPBENZ = 120.194;
Tc_NPBENZ = 638.32;
Pc_NPBENZ = 32.0;
Vc_NPBENZ = .440;
om_NPBENZ = .344;
Tb_NPBENZ = 432.2;
Ant_NPBENZ = [7.18167 1655.21 225.615]';
NPBENZ_Prop = [M_NPBENZ Tc_NPBENZ Pc_NPBENZ Vc_NPBENZ om_NPBENZ Tb_NPBENZ]';
M_NC12H26 = 170.338;
Tc_NC12H26 = 658;
Pc_NC12H26 = 18.2;
Vc_NC12H26 = .718;
om_NC12H26 = .577;
Tb_NC12H26 = 489;
Ant_NC12H26 = [7.22883 1807.47 199.381]';
NC12H26_Prop = [M_NC12H26 Tc_NC12H26 Pc_NC12H26 Vc_NC12H26 om_NC12H26 Tb_NC12H26]';
M_NC7H16 = 100.204;
Tc_NC7H16 = 540.2;
Pc_NC7H16 = 27.2;
Vc_NC7H16 = .428;
om_NC7H16 = .346;
Tb_NC7H16 = 371.5;
Ant_NC7H16 = [7.04605 1341.89 223.733]';
NC7H16_Prop = [M_NC7H16 Tc_NC7H16 Pc_NC7H16 Vc_NC7H16 om_NC7H16 Tb_NC7H16]';
M_NC10H22 = 142.285;
Tc_NC10H22 = 617.7;
Pc_NC10H22 = 20.9;
Vc_NC10H22 = .601;
om_NC10H22 = .488;
Tb_NC10H22 = 447.2;
Ant_NC10H22 = [7.21745 1693.93 216.459]';
NC10H22_Prop = [M_NC10H22 Tc_NC10H22 Pc_NC10H22 Vc_NC10H22 om_NC10H22 Tb_NC10H22]';
M_C2H5OH = 46.069;
Tc_C2H5OH = 513.92;
Pc_C2H5OH = 61.2;
Vc_C2H5OH = .168;
om_C2H5OH = .643;
Tb_C2H5OH = 351.5;
Ant_C2H5OH = [8.13484 1662.48 238.131]';
C2H5OH_Prop = [M_C2H5OH Tc_C2H5OH Pc_C2H5OH Vc_C2H5OH om_C2H5OH Tb_C2H5OH]';
M_CH3OH = 32.042;
Tc_CH3OH = 512.64;
Pc_CH3OH = 81.4;
Vc_CH3OH = .117;
om_CH3OH = .566;
Tb_CH3OH = 337.8;
Ant_CH3OH = [8.09126 1582.91 239.096]';
CH3OH_Prop = [M_CH3OH Tc_CH3OH Pc_CH3OH Vc_CH3OH om_CH3OH Tb_CH3OH]';
M_NC16H34 = 226.446;
Tc_NC16H34 = 723;
Pc_NC16H34 = 14.1;
Vc_NC16H34 = .943;
om_NC16H34 = .721;
Tb_NC16H34 = 554;
Ant_NC16H34 = [7.36235 2094.08 180.407]';
NC16H34_Prop = [M_NC16H34 Tc_NC16H34 Pc_NC16H34 Vc_NC16H34 om_NC16H34 Tb_NC16H34]';
M_H2O = 18.015;
Tc_H2O = 647.13;
Pc_H2O = 219.4;
Vc_H2O = .056;
om_H2O = .343;
Tb_H2O = 373.17;
Ant_H2O = [8.05573 1723.64 233.076]';
H2O_Prop = [M_H2O Tc_H2O Pc_H2O Vc_H2O om_H2O Tb_H2O]';

%% Data
%%
load('Data/hept50hex50_2bar.mat')
W_Comp = [W_NC16H34 W_NC7H16];
Comp_Prop = [NC16H34_Prop NC7H16_Prop];
P0 = 2;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%%
load('Data/hept50hex50_5bar.mat')
W_Comp = [W_NC16H34 W_NC7H16];
Comp_Prop = [NC16H34_Prop NC7H16_Prop];
P0 = 5;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%r=0.254
%%
load('Data/wat10hex90.mat')
W_Comp = [W_H2O W_NC16H34];
Comp_Prop = [H2O_Prop NC16H34_Prop];
P0 = 6;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)
%%
load('Data/eth50dec50_1bar.mat')
W_Comp = [W_C2H5OH W_NC10H22];
Comp_Prop = [C2H5OH_Prop NC10H22_Prop];
P0 = 1;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)
%%
load('Data/eth50dec50_10bar.mat')
W_Comp = [W_C2H5OH W_NC10H22];
Comp_Prop = [C2H5OH_Prop NC10H22_Prop];
P0 = 10;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)
%%
load('Data/meth50dod50.mat')
W_Comp = [W_NC12H26 W_CH3OH];
Comp_Prop = [NC12H26_Prop CH3OH_Prop];
P0 = 40;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%%
load('Data/meth50hex50.mat')
W_Comp = [W_NC16H34 W_CH3OH];
Comp_Prop = [NC16H34_Prop CH3OH_Prop];
P0 = 40;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%%
load('Data/meth40hex60.mat')
W_Comp = [W_NC16H34 W_CH3OH];
Comp_Prop = [NC16H34_Prop CH3OH_Prop];
P0 = 40;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%%
load('Data/hept50hex50_4bar.mat')
W_Comp = [W_NC16H34 W_NC7H16];
Comp_Prop = [NC16H34_Prop NC7H16_Prop];
P0 = 4;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%r=0.21
%%
load('Data/hept50hex50_3bar.mat')
W_Comp = [W_NC16H34 W_NC7H16];
Comp_Prop = [NC16H34_Prop NC7H16_Prop];
P0 = 3;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 0)
%%
load('Data/eth50dod50_1bar.mat')
W_Comp = [W_C2H5OH W_NC12H26];
Comp_Prop = [C2H5OH_Prop NC12H26_Prop];
P0 = 1;
SuperheatPlot(Comp_Prop, W_Comp, t, r, T, P0, 1)

%% Kinetic Model
load('Data/hept50hex50_5bar.mat')
W_Comp = [W_NC16H34 W_NC7H16];
Comp_Prop = [NC16H34_Prop NC7H16_Prop];
Antoine = [Ant_NC16H34 Ant_NC7H16];
P0 = 5;
KineticModel(Comp_Prop, Antoine, W_Comp, t, r, T, P0, 0)

%%
load('Data/eth50dec50_1bar.mat')
W_Comp = [W_C2H5OH W_NC10H22];
Comp_Prop = [C2H5OH_Prop NC10H22_Prop];
Antoine = [Ant_C2H5OH Ant_NC10H22];
P0 = 1;
KineticModel(Comp_Prop, Antoine, W_Comp, t, r, T, P0, 0)

