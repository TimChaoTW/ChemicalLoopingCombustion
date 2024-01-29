clear
clc

global L dz ro_bed porosity_bed d_p nt dti  porosity_catalyst Tf_in Pt_in Cps_t u_s  CM_CH4_in CM_H2_in CM_CO_in CM_CO2_in  CM_H2O_in Tw_in_2 CM_N2_in Ts_in
global r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 FM_t keq4 keq5 keq6 keq7 keq8 keq9 k1 k2 k3 k4 k5 k6 k7 k8 k9 ro_solid u_s_in Mwt_in Mw_t X G hf ro_t Re
global CMs_CH4_in CMs_H2_in CMs_CO_in CMs_CO2_in CMs_N2_in CMs_H2O_in CM_Ni_in CM_NiO_in V_cat CH4_kg_in 
tic

L=11; %Length of reactor, m
dz = 0.5; %Segment length

%Reactor parameters
ro_solid= 3960;%catalyst density (kg_catalyst/m^3)
ro_bed=995; %Packing density (kg_catalyst/m^3) bed, (1-porosity)*ro_cat
porosity_bed=0.4; %porosity of the bed
porosity_catalyst = 0.37; %porosity of all of the studied oxygen carriers
d_p=0.004; % diameter ofoxygen carrier[m]
nt=1; %No. of tubes
dti= 5.5; %Diameter of tube (m)
dt = 10; %Segment time
TotalTime = 1*10^-6;
Ac=pi*(dti/2)^2;
R_bar =83.151398*10^-6;
Tf_in = 723.15; %K
Pt_in = 17; %bar
%methane-fuel
% �y�t = ������n�y�v / �I���n --> us = V_dot x Ac; Ac = pi x (5.5m/2)^2 = 23.75 ���褽�� 
% ��n�y�v = ��q�y�v / �K�� --> V_dot = W_dot / ro_mass
% PV = nRT ; P:���O, V:��n, n:���ռ�, R:�z�Q����`�� R_bar = 0.08314 ��� bar x m^3 / (mol x K); T: Tf �y��ū�
% ���լy�v = �I���n  x �y�t x �@�� --> F_dot = Ac x dz x us x CM
% --> CM = F_dot / (Ac  x us)
% �K�� = ���O x ���l�q / (�z�Q����`�� x �y��?��) --> ro_mass = Pt x Mw / (R_bar x Tf); 
CH4_kg_in = 15; % �C��i�� 15���窺 CH4, ��쬰 kg / sec
CH4_Fdot_in = 15 / (16*10^-3); % CH4 �����լy�v mol/sec  --> 937.5 mol / sec
CH4_Mw = 16 *10^-3; % kg/mol
CH4_ro =  Pt_in * CH4_Mw /(R_bar*Tf_in);
CH4_Vdot= CH4_kg_in/ CH4_ro;
us_in = CH4_Vdot/Ac;

% CH4 ���@��n (CM) �p�U:
CM_CH4_in = CH4_Fdot_in/(Ac*us_in*porosity_bed); %mol/m^3 ;ro_ch4:0.72 kg/m^3,900sec

CM_H2_in = 1e-6; %mol/m^3
CM_CO_in = 1e-6; %mol/m^3
CM_CO2_in = 1e-6; %mol/m^3
CM_N2_in = 1e-6; %mol/m^3
CM_H2O_in = 1e-6; %mol/m^3
Ts_in = 925.15;%K
CMs_CH4_in = 1e-6; %mol/m^3
CMs_H2_in = 1e-6; %mol/m^3
CMs_CO_in = 1e-6; %mol/m^3
CMs_CO2_in = 1e-6; %mol/m^3
CMs_N2_in = 1e-6; %mol/m^3
CMs_H2O_in = 1e-6; %mol/m^3
CM_Ni_in = 1e-6; %kgNi/kgOC
CM_NiO_in =0.6; %kgNi/kgOC
V_cat=Ac*dz*(1-porosity_bed)
%Dynamic of heterogeneous reaction
objfcn_Dynamic = Dynamic_PBR(dt,TotalTime);
toc






