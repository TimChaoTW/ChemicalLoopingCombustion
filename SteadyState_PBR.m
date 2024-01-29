function der_c=SteadyState_PBR(t,var)

global nt dti ro_bed porosity_bed d_p porosity_catalyst L n dz ro_solid xyz 
global FM_in  Tf_in  Pt_in_2  CM CM_NiO_in Tw_read CM_Ni
global r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 k_mass h_heat U u_s Rep Re
global ro_t_in u_s_in ro_t a_specific R_bar Hr_1 Hr_2 Hr_3 Hr_4 Hr_5 Hr_6 Hr_7 Hr_8 Hr_9  Cp_t ro_t_mass a_0
global Pt_read   P_t time FM_O2_in  dtt Cps_t Tw_in_2 T_w  CM_t hf  CMt_in Mwt_in Re_p G V_cat CH4_kg_in 

global Rbar_Ts R_bar TsK TfK; % for Kgpv & Kgptc & Kgph & HR
global ReK GK porosity_bed ro_t_massK; % for Kmtc function 
global Tdiff Pdiff; % for DC.m
global Kmass; % for Kmtc.m

% for debug only, xyz is a counter for SteadyState_PBR
 xyz = xyz + 1
if xyz >=  282267
     abcde = 0;
else
     abcde = 0;
end

A_tube = nt*pi*dti^2/4; %Tube cross area (m^2)
a_specific = 6*(1-porosity_bed)/d_p; %Catalyst specific area,比表面積 by ch4 固定床反應器 (1/m);fixed bed ch6 page 11 
R_atm = 0.08206; % Gas constant (atm*m^3/(kmol*K))
R_bar = 8.3151398; % Gas constant (Pa*m^3/(mol*K))
CM_NiO_in = 0.6; % initial concentration
CM_Ni_in = 1e-6; % initial concentration
a_0 = 10000;%specific surface area of the oxygen carrier(m^2/kg)
FM_in = 937.5;%mol/s
CMt_in = 282.72;%mol/m^3
u_s_in = FM_in/(CMt_in*A_tube); %Gas superficial velocity (m/s)

Mw_CH4=16/1000; %kg/mol
Mw_H2=2/1000; %kg/mol
Mw_CO=28/1000; %kg/mol
Mw_CO2=44/1000; %kg/mol
Mw_H2O=18/1000; %kg/mol
Mw_O2=32/1000; %kg/mol
Mw_N2=28/1000;%kg/mol

i= 1:n;
CM(1,i) = var(1:n); %CH4 mol/m^3
CM(2,i) = var(n+1:2*n); %H2 mol/m^3
CM(3,i) = var(2*n+1:3*n); %CO mol/m^3
CM(4,i) = var(3*n+1:4*n); %CO2 mol/m^3
CM(5,i) = var(4*n+1:5*n); %N2 mol/m^3
CM(6,i) = var(5*n+1:6*n); %H2O mol/m^3
Tf(i) = var(6*n+1:7*n);
Ts(i) = var(7*n+1:8*n); 
CMs(1,i) = var(8*n+1:9*n); %CH4
CMs(2,i) = var(9*n+1:10*n); %H2
CMs(3,i) = var(10*n+1:11*n); %CO
CMs(4,i) = var(11*n+1:12*n); %CO2
CMs(5,i) = var(12*n+1:13*n); %N2
CMs(6,i) = var(13*n+1:14*n); %H2O
CM_Ni(i) = var(14*n+1:15*n); %concentration of Ni kgNi/kgOC
CM_NiO(i) = var(15*n+1:16*n); %concentration of NiO kgNi/kgOC
X(i) = var(16*n+1:17*n); %conversion of oxygen carrier
 
 
 for i=1:n
CM_t(i) = CM(1,i)+CM(2,i)+CM(3,i)+CM(4,i)+CM(5,i)+CM(6,i); %mol/m^3    
 end
 
for i=1:n
Pt(i) =17-(i-1)*0.1;    
end
 
%Molecular weight
Mw_CH4=16/1000; % Molecular weight of CH4 kg/mol 重量分子量
Mw_H2=2/1000;   % Molecular weight of H2 kg/mol 重量分子量
Mw_CO=28/1000;  % Molecular weight of CO kg/mol 重量分子量
Mw_CO2=44/1000; % Molecular weight of CO2 kg/mol 重量分子量
Mw_O2=32/1000;  % Molecular weight of O2 kg/mol 重量分子量
Mw_H2O=18/1000; % Molecular weight of H2O kg/mol 重量分子量
Mw_NiO=74.69/1000;% Molecular weight of NiO kg/mol 重量分子量
Mw_N2=28/1000; %Molecular weight of N2 kg/mol 重量分子量
Mw_Ni=58.69/1000; %Molecular weight of Ni kg/mol 重量分子量

for i=1:n
%Mole fraction
y(1,i) = CM(1,i)/CM_t(i); %CH4 莫耳濃度比
y(2,i) = CM(2,i)/CM_t(i); %H2 莫耳濃度比
y(3,i) = CM(3,i)/CM_t(i); %CO 莫耳濃度比
y(4,i) = CM(4,i)/CM_t(i); %CO2 莫耳濃度比
y(5,i) = CM(5,i)/CM_t(i); %N2 莫耳濃度比
y(6,i) = CM(6,i)/CM_t(i); %H2O 莫耳濃度比
Mw_t(i) = y(1,i)*Mw_CH4 + y(2,i)*Mw_H2 + y(3,i)*Mw_CO + y(4,i)*Mw_CO2 + y(5,i)*Mw_N2 + y(6,i)*Mw_H2O; % kg/mol 總分子量
ro_t_mass(i) = CM_t(i)*Mw_t(i); %Vapor density(kg/m^3)
end

for i= 1:n
% TfK = Tf(i);
    %Fluid vicosity form line 94 to 114 but no H2O,Unit kg/(m*s), viscosity form  Paper Heterogenous page 5 table 4
if Tf(i) > 1000
    mu_CH4(i) = (5.255e-7*1000^0.59006/(1+105.67/1000)) ; %kg/(m*s) CH4 黏度
else
    mu_CH4(i) = (5.255e-7*Tf(i)^0.59006/(1+105.67/Tf(i))); %kg/(m*s) CH4 黏度　
end

if Tf(i) > 1250
    mu_CO(i) = (1.113e-6*1250^0.5338/(1+94.7/1250)); %kg/(m*s) CO 黏度
else
    mu_CO(i) = (1.113e-6*Tf(i)^0.5338/(1+94.7/Tf(i))); %kg/(m*s) CO 黏度
end

if Tf(i) > 1073
    mu_H2O(i) = (1.71e-8*1073^1.1146); %kg/(m*s) H2O 黏度
else
    mu_H2O(i) = (1.71e-8*Tf(i)^1.1146); %kg/(m*s) H2O 黏度
end

mu_H2(i) = (1.797e-7*Tf(i)^0.685/(1+(-0.59)/Tf(i)+140/Tf(i)^2)); %kg/(m*s) H2 黏度
mu_CO2(i) = Kgptc(2.148e-6,0.46,290,0,Tf(i)); %kg/(m*s)  CO2 黏度
mu_N2(i) = Kgptc(6.5592e-7,0.6081,54.714,0,Tf(i)); %kg/(m*s) N2 黏度
mu_O2(i) = Kgptc(1.101e-6,0.5634,96.3,0,Tf(i)); %kg/(m*s) O2 黏度
mu_t(i) = y(1,i)*mu_CH4(i)+y(2,i)*mu_H2(i)+y(3,i)*mu_CO(i)+y(4,i)*mu_CO2(i)+y(5,i)*mu_N2(i)+y(6,i)*mu_H2O(i); %kg/(m*s) 總氣黏度

%------------------------------------------------------------------------------------------------
%Fluid heat capacity; %line 119 to 126 from dynamic page 256 Table 12.3; check OK 2018/11/24
% Cp_CH4_old = 33298 + 79933*(2086.9/(Tf(i)*sinh((2086.9/Tf(i))^2))) + 41602*(991.96/(Tf(i)*sinh((991.96/Tf(i))^2)))*1e-3; %J/(mol*K)
% Cp_H2_old = 27617 + 9560*(2466/(Tf(i)*sinh((2466/Tf(i))^2))) + 3760*(567.6/(Tf(i)*sinh((567.6/Tf(i))^2)))*1e-3 ;%J/(mol*K)
% Cp_H2O_old = 33363 + 26790*(2610.5/(Tf(i)*sinh((2610.5/Tf(i))^2))) + 8896*(1169/(Tf(i)*sinh((1169/Tf(i))^2)))*1e-3; %J/(mol*K)
% Cp_CO_old = 29108 + 8773*(3085.1/(Tf(i)*sinh((3085.1/Tf(i))^2))) + 8455.3*(1538.2/(Tf(i)*sinh((1538.2/Tf(i))^2)))*1e-3; %J/(mol*K)
% Cp_CO2_old = 29370 + 34540*(1428/(Tf(i)*sinh((1428/Tf(i))^2))) + 26400*(588/(Tf(i)*sinh((588/Tf(i))^2)))*1e-3; %J/(mol*K)
% Cp_N2_old = 29105 + 8614.9*(1701.6/(Tf(i)*sinh((1701.6/Tf(i))^2))) + 103.47*(909.79/(Tf(i)*sinh((909.79/Tf(i))^2)))*1e-3;% J/mol/K
% Cp_O2_old = 29103+10040*(2526.5/(Tf(i)*sinh((2526.5/Tf(i))^2))) + 9356*(1153.8/(Tf(i)*sinh((1153.8/Tf(i))^2)))*1e-3 ;%J/(mol*K)
% Cpi = A + B*(C/(T*sinh((C/T)^2))) + D*(E/(T*sinh((E/T)^2)))*1e-3; T = TfK = Tf(i)
Cp_CH4(i) = Kgphc(33298,79933,2086.9,41602,991.96,Tf(i)); %J/(mol*K) CH4 莫耳熱容量 (氣體比熱)
Cp_H2(i) = Kgphc(27617,9560,2466,3760,567.6,Tf(i)); %J/(mol*K) H2 莫耳熱容量 (氣體比熱)
Cp_H2O(i) = Kgphc(33363,26790,2610.5,8896,1169,Tf(i)); %J/(mol*K) H2O 莫耳熱容量 (氣體比熱)
Cp_CO(i) = Kgphc(29108,8773,3085.1,8455.3,1538.2,Tf(i)); %J/(mol*K) CO 莫耳熱容量 (氣體比熱)
Cp_CO2(i) = Kgphc(29370,34540,1428,26400,588,Tf(i)); %J/(mol*K) CO2 莫耳熱容量 (氣體比熱)
Cp_N2(i) = Kgphc(29105,8614.9,1701.6,103.47,909.79,Tf(i));% J/mol/K N2 莫耳熱容量 (氣體比熱)
Cp_O2(i) = Kgphc(29103,10040,2526.5,9356,1153.8,Tf(i)); %J/(mol*K) O2 莫耳熱容量 (氣體比熱)
Cp_t(i) = (y(1,i)*Cp_CH4(i)+y(2,i)*Cp_H2(i)+y(3,i)*Cp_CO(i)+y(4,i)*Cp_CO2(i)+y(5,i)*Cp_N2(i)+y(6,i)*Cp_H2O(i)); %J/(mol*K) 總莫耳熱容量 (氣體比熱)

%------------------------------------------------------------------------------------------------
%Solid heat capacity, >=-> See Dynamic page 259 Table 12.6's Line Nio's Cp0 Cp1 Cp2; check OK 2018/11/24
CpsA = 790*Mw_NiO;
CpsB = -2.06e-1*Mw_NiO;
CpsC = 1.43e-4*Mw_NiO;
Cps_t(i) = CpsA + CpsB*Tf(i) + CpsC*Tf(i)^2; %J/(mol*K) NiO 觸煤比熱，或是固體比熱

%------------------------------------------------------------------------------------------------
%Reduction from paper supplementary page 13 table SI.2 line 137 to 145
% r1: H2 + NiO --> Ni + H2O
% r2: CO + NiO --> Ni + CO2
% r3: CH4 + NiO --> Ni + 2 H2 + CO
% r4: CH4 + H2O <--> CO + 3 H2
% r5: CO + H2O <--> CO2 + H2
% r6: CH4 + CO2 <--> 2 CO + 2 H2
% r7: CH4  <--> C + 2 H2
% r8: C +  H2O <--> CO + H2
% r9: C + CO2 <--> 2 CO

%------------------------------------------------------------------------------------------------
%Heterogeneous equilibrium constant & adsorption coefficient, line 149 to 165 ,Paper form model page10-table 5 ; check OK 2018/11/24
% keq4_old = 1.2e13*exp(-223/(R_bar*Ts(i))); % R6 bar^2 for r4
% keq5_old = 1.77e-02*exp(36.6/(R_bar*Ts(i))); % R7 for r5
% keq6_old = 6.78e14*exp(-260/(R_bar*Ts(i))); % R8 bar^2 for r6
% keq7_old = 2.98e5*exp(-84.4/(R_bar*Ts(i))); % R9 bar^2 for r7
% keq8_old = 4.02e7*exp(-139/(R_bar*Ts(i))); % R10 bar^1 for r8
% keq9_old = 2.28e9*exp(-175/(R_bar*Ts(i))); % R11 bar^1 for r9

% Kj = Kj0*exp(-1*DeltaH/(R_bar*TsK))　
Rbar_Ts = R_bar * Ts(i); % 氣體常數 R 的單位則使用 8.314 J·mol-1·K-1, 出自 '化學平衡常數有?有單位？.pdf' 第 6 頁
keq4(i) = 1.2e13*exp(-223/Rbar_Ts); % R6 [r4: CH4 + H2O <--> CO + 3 H2 平衡常數]
keq5(i) = 1.77e-02*exp(36.6/Rbar_Ts); % R7 [r5: CO + H2O <--> CO2 + H2 平衡常數]
keq6(i) = 6.78e14*exp(-260/Rbar_Ts); % R8 [r6: CH4 + CO2 <--> 2 CO + 2 H2 平衡常數]
keq7(i) = 2.98e5*exp(-84.4/Rbar_Ts); % R9 [r7: CH4  <--> C + 2 H2 平衡常數]
keq8(i) = 4.02e7*exp(-139/Rbar_Ts); % R10 [r8: C +  H2O <--> CO + H2 平衡常數]
keq9(i) = 2.28e9*exp(-175/Rbar_Ts); % R11 [r9: C + CO2 <--> 2 CO 平衡常數]

%------------------------------------------------------------------------------------------------
%Heterogeneous adsorption parameters of component, line 160 to 173 from supplementary page 18 table SI.8; check OK 2018/11/24
% kch4_rd45_old = 6.65e-4*exp(38.3/R_bar/Ts(i)); % Kch4 R4 & R5 bar^-1
% kco_rd45_old = 8.23e-5*exp(70.6/R_bar/Ts(i)); % Kco R4 & R5 bar^-1
% kh2_rd45_old = 6.12e-9*exp(82.9/R_bar/Ts(i)); %bar^-1
% kh2o_rd45_old = 1.77e5*exp(-88.7/R_bar/Ts(i)); %bar^-1
% kch4_rd6_old = 4.04e-4*exp(74.6/R_bar/Ts(i)); %bar^-1
% kch4_rd7_old = 2.1e-1*exp(-0.567/R_bar/Ts(i)); %bar^-1
% kh2_rd7_old = 5.18e7*exp(-133/R_bar/Ts(i)); %bar^-1
% kch4_rd8_old = 3.49; %bar^-1
% kco_rd8_old = 4.73e-6*exp(97.8/R_bar/Ts(i)) ;%bar^-1
% kh2_rd8_old = 1.83e13*exp(-216/R_bar/Ts(i)) ;%bar^-1
% kco_rd9_old = 7.34e-6*exp(100.4/R_bar/Ts(i)) ;%bar^-1
% kco2_rd9_old = 8.17e7*exp(-104/R_bar/Ts(i)) ;%bar
% kco_rd13_old = 5.22*exp(-8.37/R_bar/Ts(i)) ;%bar^-1
% Kij =  Kij0*exp(-DeltaHij0/(R*T))

kch4_rd45(i) = 6.65e-4*exp(38.3/Rbar_Ts); % Kch4 R4 & R5 bar^-1　吸附常數　
kco_rd45(i) = 8.23e-5*exp(70.6/Rbar_Ts); % Kco R4 & R5 bar^-1
kh2_rd45(i) = 6.12e-9*exp(82.9/Rbar_Ts); % Kh2 R4 & R5 bar^-1
kh2o_rd45(i) = 1.77e5*exp(-88.7/Rbar_Ts); % Kh2o R4 & R5 bar^-1
kch4_rd6(i) = 4.04e-4*exp(74.6/Rbar_Ts); % Kch4 R6 bar^-1
kch4_rd7(i) = 2.1e-1*exp(-0.567/Rbar_Ts); % Kch4 R7 bar^-1
kh2_rd7(i) = 5.18e7*exp(-133/Rbar_Ts); % Kh2 R7 bar^-1
kch4_rd8(i) = 3.49; % Kch4 R8 bar^-1
kco_rd8(i) = 4.73e-6*exp(97.8/Rbar_Ts) ;% Kco R8 bar^-1
kh2_rd8(i) = 1.83e13*exp(-216/Rbar_Ts) ;% Kh2 R8 bar^-1
kco_rd9(i) = 7.34e-6*exp(100.4/Rbar_Ts) ;% Kco R9 bar^-1
kco2_rd9(i) = 8.17e7*exp(-104/Rbar_Ts) ;% Kco2 R9 bar
kco_rd13(i) = 5.22*exp(-8.37/Rbar_Ts) ;% Kco R13 bar^-1

%-----------------------------------------------------------------------------------------------
%Heterogeneous reaction rates 
%ｂy pass the real value using heterogenous table 10 page 8
%line 177 to 185 from dynamic page 40equation 2.2 & supplementary page 17 table SI.5, check OK 2018/11/24
% k1_old = 7.7e-2*exp(-27/R_bar/Ts(i)); %m/s
% k2_old = 4.22e-4*exp(-37/R_bar/Ts(i)); %m/s
% k3_old = 2.32e-4*exp(-39/R_bar/Ts(i)); %m/s
% k4_old = 1.01e5*exp(-168/R_bar/Ts(i)) ;%bar^0.5*kmol/kgNi*s
% k5_old = 7.25e1*exp(-102/R_bar/Ts(i)) ;%mol/bar*kgNi*s
% k6_old = 9.21e-1*exp(-124/R_bar/Ts(i)); %mol/bar*kgNi*s
% k7_old = 7e-1*exp(-59/R_bar/Ts(i)); %bar^-1
% k8_old = 7.1e-1*exp(-132/R_bar/Ts(i)); %mol/kgNi*s
% k9_old = 1.47e13*exp(-365/R_bar/Ts(i)); %mol/kgNi*s
% Kj = Kj0*exp(-1*DeltaH/(R*T))
% k1 = 7.7e-2*exp(-27/Rbar_Ts); % R1 m/s
% k2 = 4.22e-4*exp(-37/Rbar_Ts); % R2 m/s
% k3 = 2.32e-4*exp(-39/Rbar_Ts); % R3 m/s
% k4 = 1.01e5*exp(-168/Rbar_Ts) ;% R4 bar^0.5*kmol/kgNi*s
% k5 = 7.25e1*exp(-102/Rbar_Ts) ;% R5 mol/bar*kgNi*s
% k6 = 9.21e-1*exp(-124/Rbar_Ts); % R6 mol/bar*kgNi*s
% k7 = 7e-1*exp(-59/Rbar_Ts); % R7 bar^-1
% k8 = 7.1e-1*exp(-132/Rbar_Ts); % R8 mol/kgNi*s
% k9 = 1.47e13*exp(-365/Rbar_Ts); % R9 mol/kgNi*s
k1=2.26e-05;
k2=1.32e-05;
k3=8.16e-06;
k4=5.44e-03;
k5=1.68e-05;
k6=9.37e-05;
k7=6.07e-05;
k8=5.40e-07;
k9=1.07e-04;

%-----------------------------------------------------------------------------------------------
%reduction phase
deno1(i) = (1+kco_rd45(i)*Pt(i)*y(3,i)+kh2_rd45(i)*Pt(i)*y(2,i)+kch4_rd45(i)*Pt(i)*y(1,i)+(kh2o_rd45(i)*Pt(i)*y(6,i)/(Pt(i)*y(2,i))))^2;

%for simplify the equation r1-r3,paper for supplementary page 16 Table SI.4 
% r1(i) = a_0*k1/(Pt(i)^1.39)*CM(2,i)*CM_NiO_in*(1-X(i))*(-log(1-X(i))); %mol/kgcat-s
% r2(i) = a_0*k2/(Pt(i)^1.21)*CM(3,i)*CM_NiO_in*(1-X(i))*(-log(1-X(i)));%mol/kgcat-s
% r3(i) = a_0*k3/(Pt(i)^1.01)*CM(1,i)*CM_NiO_in*(1-X(i))*(-log(1-X(i))); %mol/kgcat-s


r1(i) = a_0*k1*CM(2,i)*CM_NiO_in*(1-X(i));%mol/kgcat-s
r2(i) = a_0*k2*CM(3,i)*CM_NiO_in*(1-X(i)); %mol/kgcat-s
r3(i) = a_0*k3*CM(1,i)*CM_NiO_in*(1-X(i))*CM_Ni(i); %mol/kgcat-s

r4(i) = (k4/(Pt(i)*y(2,i))^2.5)*(Pt(i)*y(1,i)*Pt(i)*y(6,i)-((Pt(i)*y(2,i))^3*Pt(i)*y(3,i)/keq4(i)))/deno1(i); %mol/kgOC-s
r5(i) = (k5/(Pt(i)*y(2,i)))*(Pt(i)*y(3,i)*Pt(i)*y(6,i)-(Pt(i)*y(2,i)*Pt(i)*y(4,i)/keq5(i)))/deno1(i); %mol/kgOC-s
r6(i) = k6*kch4_rd6(i)*(Pt(i)*y(1,i)*Pt(i)*y(4,i)-((Pt(i)*y(2,i))^2*(Pt(i)*y(3,i))^2/keq6(i)))/(1+kch4_rd6(i)*Pt(i)*y(1,i)); %mol/kgOC-s
r7(i) = k7*kch4_rd7(i)*(Pt(i)*y(1,i)-((Pt(i)*y(2,i))^2/keq7(i)))/(1+((Pt(i)*y(2,i))^1.5/kch4_rd7(i))+kch4_rd7(i)*Pt(i)*y(1,i))^2 ;%mol/kgOC-s
r8(i) = k8*(Pt(i)*y(6,i)/(Pt(i)*y(2,i))-(Pt(i)*y(3,i)/keq8(i)))/(1+(kch4_rd8(i)*Pt(i)*y(1,i))+(Pt(i)*y(6,i)/Pt(i)*y(2,i))+((Pt(i)*y(2,i))^1.5/kh2_rd8(i)))^2 ;%mol/kgOC-s
r9(i) = k9/(kco_rd9(i)*kco2_rd9(i))*((Pt(i)*y(4,i)/(Pt(i)*y(3,i)))-(Pt(i)*y(3,i)/keq9(i)))/(1+(Pt(i)*y(4,i)/kco_rd9(i)/kco2_rd9(i)/Pt(i)*y(3,i))+(kco_rd9(i)*Pt(i)*y(3,i)))^2; %mol/kgOC-s

end

%Diffusion coefficients of vapor(m^2/s) 
%form Paper Heterogenous page 4 table 1 ???
M_CH4_H2O = 2*(1/Mw_CH4+1/Mw_H2O)^(-1); %2[(1 /MA)+(1 /MB)]^-1
M_CH4_H2 = 2*(1/Mw_CH4+1/Mw_H2)^(-1);%2[(1 /MA)+(1 /MB)]^-1
M_CH4_CO2 = 2*(1/Mw_CH4+1/Mw_CO2)^(-1); %2[(1 /MA)+(1 /MB)]^-1
M_CH4_CO = 2*(1/Mw_CH4+1/Mw_CO)^(-1) ;%2[(1 /MA)+(1 /MB)]^-1
M_CH4_N2 = 2*(1/Mw_CH4+1/Mw_N2)^(-1) ;%2[(1 /MA)+(1 /MB)]^-1

M_H2O_H2 = 2*(1/Mw_H2O+1/Mw_H2)^(-1);%2[(1 /MA)+(1 /MB)]^-1
M_H2O_CO2 = 2*(1/Mw_H2O+1/Mw_CO2)^(-1); %2[(1 /MA)+(1 /MB)]^-1
M_H2O_CO = 2*(1/Mw_H2O+1/Mw_CO)^(-1) ;%2[(1 /MA)+(1 /MB)]^-1
M_H2O_N2 = 2*(1/Mw_H2O+1/Mw_N2)^(-1); %2[(1 /MA)+(1 /MB)]^-1

M_H2_CO2 = 2*(1/Mw_H2+1/Mw_CO2)^(-1) ;%2[(1 /MA)+(1 /MB)]^-1
M_H2_CO = 2*(1/Mw_H2+1/Mw_CO)^(-1) ;%2[(1 /MA)+(1 /MB)]^-1
M_H2_N2 = 2*(1/Mw_H2+1/Mw_N2)^(-1) ;%2[(1 /MA)+(1 /MB)]^-1

M_CO2_CO = 2*(1/Mw_CO2+1/Mw_CO)^(-1); %2[(1 /MA)+(1 /MB)]^-1
M_CO2_N2 = 2*(1/Mw_CO2+1/Mw_N2)^(-1); %2[(1 /MA)+(1 /MB)]^-1

M_CO_N2 = 2*(1/Mw_CO+1/Mw_N2)^(-1); %2[(1 /MA)+(1 /MB)]^-1

sigma_CH4 = 25.14; %Diffusion Volumes of Simple Molecules
sigma_H2O = 13.1; %Diffusion Volumes of Simple Molecules
sigma_H2 = 6.12; %Diffusion Volumes of Simple Molecules
sigma_CO2 = 26.9; %Diffusion Volumes of Simple Molecules
sigma_CO = 18; %Diffusion Volumes of Simple Molecules
sigma_O2 = 16.3; %Diffusion Volumes of Simple Molecules
sigma_N2 = 17.9; %Diffusion Volumes of Simple Molecules

%initial superficial velocity
u_s(1) = u_s_in;%m/s
CM_t(1) = CMt_in ;%mol/m^3

for i = 2:n
%superficial velocity
u_s(i) = CH4_kg_in /(Mw_t(i)*CM_t(i)*A_tube);%Gas superficial velocity (m/s)

end

for i=1:n
%-------------------------------------------------------------------------------------------------------------------------------------
Tdiff = 0.00143*Tf(i)^1.75*1e-4;
Pdiff = Pt(i)/101330*1.013;
D_CH4_H2O(i) = DC(M_CH4_H2O,sigma_CH4,sigma_H2O); %DIFFUSION COEFFICIENTS(m^2/s)
D_CH4_H2(i) = DC(M_CH4_H2,sigma_CH4,sigma_H2); %DIFFUSION COEFFICIENTS(m^2/s)
D_CH4_CO2(i) = DC(M_CH4_CO2,sigma_CH4,sigma_CO2); %DIFFUSION COEFFICIENTS(m^2/s)
D_CH4_CO(i) = DC(M_CH4_CO,sigma_CH4,sigma_CO); %DIFFUSION COEFFICIENTS(m^2/s)
D_CH4_N2(i) = DC(M_CH4_N2,sigma_CH4,sigma_N2); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2O_H2(i) = DC(M_H2O_H2,sigma_H2O,sigma_H2); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2O_CO2(i) = DC(M_H2O_CO2,sigma_H2O,sigma_CO2); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2O_CO(i) = DC(M_H2O_CO,sigma_H2O,sigma_CO); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2O_N2(i) = DC(M_H2O_N2,sigma_H2O,sigma_N2); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2_CO2(i) = DC(M_H2_CO2,sigma_H2,sigma_CO2); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2_CO(i) = DC(M_H2_CO,sigma_H2,sigma_CO); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2_N2(i) = DC(M_H2_N2,sigma_H2,sigma_N2); %DIFFUSION COEFFICIENTS(m^2/s)
D_CO2_CO(i) = DC(M_CO2_CO,sigma_CO2,sigma_CO); %DIFFUSION COEFFICIENTS(m^2/s)
D_CO2_N2(i) = DC(M_CO2_N2,sigma_CO2,sigma_N2); %DIFFUSION COEFFICIENTS(m^2/s)
D_CO_N2(i) = DC(M_CO_N2,sigma_CO,sigma_N2); %DIFFUSION COEFFICIENTS(m^2/s)

D_CH4(i) = (1-y(1,i))*(y(6,i)/D_CH4_H2O(i) + y(2,i)/D_CH4_H2(i) + y(4,i)/D_CH4_CO2(i) + y(3,i)/D_CH4_CO(i)+y(5,i)/D_CH4_N2(i))^(-1); %DIFFUSION COEFFICIENTS(m^2/s)
D_H2O(i) = (1-y(6,i))*(y(1,i)/D_CH4_H2O(i) + y(2,i)/D_H2O_H2(i) + y(4,i)/D_H2O_CO2(i) + y(3,i)/D_H2O_CO(i)+y(5,i)/D_H2O_N2(i))^(-1) ;%DIFFUSION COEFFICIENTS(m^2/s)
D_H2(i) = (1-y(2,i))*(y(6,i)/D_H2O_H2(i) + y(1,i)/D_CH4_H2(i) + y(4,i)/D_H2_CO2(i) + y(3,i)/D_H2_CO(i)+y(5,i)/D_H2_N2(i))^(-1); %DIFFUSION COEFFICIENTS(m^2/s)
D_CO2(i) = (1-y(4,i))*(y(6,i)/D_H2O_CO2(i) + y(2,i)/D_H2_CO2(i) + y(1,i)/D_CH4_CO2(i) + y(3,i)/D_CO2_CO(i)+y(5,i)/D_CO2_N2(i))^(-1); %DIFFUSION COEFFICIENTS(m^2/s)
D_CO(i) = (1-y(3,i))*(y(6,i)/D_H2O_CO(i) + y(2,i)/D_H2_CO(i) + y(4,i)/D_CO2_CO(i) + y(1,i)/D_CH4_CO(i)+y(5,i)/D_CO_N2(i))^(-1) ;%DIFFUSION COEFFICIENTS(m^2/s)
D_N2(i) = (1-y(5,i))*(y(6,i)/D_H2O_N2(i) + y(2,i)/D_H2_N2(i) + y(4,i)/D_CO2_N2(i) + y(1,i)/D_CH4_N2(i)+y(3,i)/D_CO_N2(i))^(-1); %DIFFUSION COEFFICIENTS(m^2/s)

%-------------------------------------------------------------------------------------------------------------------------------------
%Thermal conductivity line 279 to 286, form Paper Heterogenous page 6 table 5 check Ok but O2 unknow 2018/11/24
% ln_CH4_old = (8.398e-6*Tf(i)^1.4268/(1+(-49.654/Tf(i)))) ;%J/(s*m*K)
% ln_H2_old = (2.653e-3*Tf(i)^0.7452/(1+(12/Tf(i)))); %J/(s*m*K)
% ln_H2O_old = (6.2041e-6*Tf(i)^1.3973); %J/(s*m*K)
% ln_CO_old = (5.9882e-4*Tf(i)^0.6863/(1+(57.13/Tf(i))+(501.92/Tf(i))^2)); %J/(s*m*K)
% ln_CO2_old = (3.69*Tf(i)^(-0.3838)/(1+(964/Tf(i))+(1.86e6/Tf(i))^2)); %J/(s*m*K)
% ln_N2_old = (3.3143e-4*Tf(i)^0.7722/(1+(16.323/Tf(i))+(373.72/Tf(i))^2)); %J/(s*m*K)
% ln_O2_old = 0.00044994*Tf(i)^0.7456/(1+56.699/Tf(i)); %KJ/(s*m*K)
% Lambda = A*T^B / ( 1 + C/T + (D/T^2));
ln_CH4(i) = Kgptc(8.398e-6,1.4268,-49.654,0,Tf(i)); %J/(s*m*K)
ln_H2(i) = Kgptc(2.653e-3,0.7452,12,0,Tf(i)); %J/(s*m*K)
ln_H2O(i) = Kgptc(6.2041e-6,1.3973,0,0,Tf(i)); %J/(s*m*K)
ln_CO(i) = Kgptc(5.9882e-4,0.6863,57.13,501.92,Tf(i)); %J/(s*m*K)
ln_CO2(i) = Kgptc(3.69,-0.3838,964,1.86e6,Tf(i)); %J/(s*m*K)
ln_N2(i) = Kgptc(3.3143e-4,0.7722,16.323,373.72,Tf(i)); %J/(s*m*K)
ln_O2(i) = Kgptc(4.4994e-4,0.7456,56.699,0,Tf(i)); %KJ/(s*m*K)
ln_t(i) = y(1,i)*ln_CH4(i)+y(2,i)*ln_H2(i)+y(3,i)*ln_CO(i)+y(4,i)*ln_CO2(i)+y(5,i)*ln_N2(i)+y(6,i)*ln_H2O(i); %J/(s*m*K)

%-------------------------------------------------------------------------------------------------------------------------------------
%Prandtl number ; Pr=cp*mu/ln from wiki,unit Cp (J/kg*K) * mu [N(kg*m/s^2)*s/m^2] / k [W(J/s)/(m*K)] check ok 2018/11/24
Pr_CH4(i) = Cp_CH4(i)*mu_CH4(i)/ln_CH4(i);
Pr_H2(i) = Cp_H2(i)*mu_H2(i)/ln_H2(i);
Pr_CO(i) = Cp_CO(i)*mu_CO(i)/ln_CO(i);
Pr_CO2(i) = Cp_CO2(i)*mu_CO2(i)/ln_CO2(i);
Pr_H2O(i) = Cp_H2O(i)*mu_H2O(i)/ln_H2O(i);
Pr_O2(i) = Cp_O2(i)*mu_O2(i)/ln_O2(i);
Pr_N2(i) = Cp_N2(i)*mu_N2(i)/ln_N2(i);
Pr_t(i) = y(1,i)*Pr_CH4(i)+y(2,i)*Pr_H2(i)+y(3,i)*Pr_CO(i)+y(4,i)*Pr_CO2(i)+y(5,i)*Pr_N2(i)+y(6,i)*Pr_H2O(i);
    
ro_t(i) = ro_t_mass(i) ;%Molae volume(mol/m^3)
Re(i) = u_s(i)*ro_t_mass(i)*d_p/mu_t(i);
Re_p(i)=ro_t_mass(i)*u_s(i)*d_p/(1-porosity_bed)/mu_t(i); %form heterogenous page 4 equation 22
end

% for i = 2:n
% %pressure hetergeneous page 4 function 21
% % Pt(i) = -(150/Re_p(i)+1.75)*(1-porosity_bed)/porosity_bed^3/d_p*ro_t_mass(i)*u_s(i)^2; %moment blance 
% Pt(i) = -(150/Re_p(i)+1.75)*(1-porosity_bed)/porosity_bed^3/d_p*ro_t_mass(i)*u_s(i)^2*dz ;%moment blance 
% end

for i = 1:n
    
%-------------------------------------------------------------------------------------------------------------------------------------
%Schmidt number ; mu / (ro * mass_diffusivity ) from wiki ==> mu (kg/m*s) / [ ro (kg/m^3) * mass_diffusivity (m^2/s) ] check ok 2018/11/24
Sc_CH4(i) = mu_CH4(i)/ro_t(i)/D_CH4(i);
Sc_H2(i) = mu_H2(i)/ro_t(i)/D_H2(i);
Sc_CO(i) = mu_CO(i)/ro_t(i)/D_CO(i);
Sc_CO2(i) = mu_CO2(i)/ro_t(i)/D_CO2(i);
Sc_H2O(i) = mu_H2O(i)/ro_t(i)/D_H2O(i);
Sc_N2(i) = mu_N2(i)/ro_t(i)/D_N2(i);

Sc_t(i) = y(1,i)*Sc_CH4(i)+y(2,i)*Sc_H2(i)+y(3,i)*Sc_CO(i)+y(4,i)*Sc_CO2(i)+y(5,i)*Sc_N2(i)+y(6,i)*Sc_H2O(i);
%-------------------------------------------------------------------------------------------------------------------------------------
%mass flux of the gas phase
% G(i) = ro_t_mass(i)*u_s(i);%kg/m^2/s
G(i) = ro_t(i)*u_s(i);%kg/m^2/s
%-------------------------------------------------------------------------------------------------------------------------------------
%heat transfercoefficient between and oxygen carrier particle; line 330 from dynamic page 59 equation 4.32; check OK 2018/11/24
hf(i) = 1.37*(0.357/porosity_bed)*Cp_t(i)/Mw_t(i)*G(i)*Re(i)^-0.359*Sc_t(i)^-(2/3); % W/m2/K

%-------------------------------------------------------------------------------------------------------------------------------------
%Mass transfer coefficient(m/s); line 331 to 336 from dynamic page 59 equation 4.31; check OK 2018/11/24
% k_mass1_old = 0.357*Re(i)^-0.359*Sc_CH4(i)^-(2/3)*(G(i)/(porosity_bed*ro_t_mass(i)));
% k_mass2_old = 0.357*Re(i)^-0.359*Sc_H2(i)^-(2/3)*(G(i)/(porosity_bed*ro_t_mass(i)));
% k_mass3_old = 0.357*Re(i)^-0.359*Sc_CO(i)^-(2/3)*(G(i)/(porosity_bed*ro_t_mass(i)));
% k_mass5_old = 0.357*Re(i)^-0.359*Sc_N2(i)^-(2/3)*(G(i)/(porosity_bed*ro_t_mass(i)));
% k_mass6_old = 0.357*Re(i)^-0.359*Sc_H2O(i)^-(2/3)*(G(i)/(porosity_bed*ro_t_mass(i)));

% GK = G(i);
% ro_t_massK = ro_t_mass(i);
% ReK = Re(i);
% k_mass(1,i) = Kmtc(Sc_CH4(i));
% k_mass(2,i) = Kmtc(Sc_H2(i));
% k_mass(3,i) = Kmtc(Sc_CO(i));
% k_mass(5,i) = Kmtc(Sc_N2(i));
% k_mass(6,i) = Kmtc(Sc_H2O(i));

Kmass = 0.357*Re(i)^-0.359*(G(i)/(porosity_bed*ro_t_mass(i)));
k_mass(1,i) = Kmass*Sc_CH4(i)^-(2/3);
k_mass(2,i) = Kmass*Sc_H2(i)^-(2/3);
k_mass(3,i) = Kmass*Sc_CO(i)^-(2/3);
k_mass(5,i) = Kmass*Sc_N2(i)^-(2/3);
k_mass(6,i) = Kmass*Sc_H2O(i)^-(2/3);

%-------------------------------------------------------------------------------------------------------------------------------------
% %Heterogeneous heat of reaction

%Heterogeneous heat of reaction; See 'H計算.xlsx'
% DeltaH = Redu*1e3 + (a0*(TsK-RT)+a1/2*(TsK^2-RT^2)+a2/3*(TsK^3-RT^3)+a3/4*(TsK^4-RT^4)+a4/5*(TsK^5-RT^5)+a5/6*(TsK^6-RT^6))*1e-3 ; %J/mol
% TsK = Ts(i);
Hr_1(i) = HR(169.91,1.41e5,2.69e1,-2.96e-1,4.34e-4,-2.53e-7,5.4e-11,Ts(i)); %J/mol
Hr_2(i) = HR(-41.5,7.06e4,-1.79e2,3.54e-1,-3.22e-4,1.47e-7,-2.65e-11,Ts(i)); %J/mol
Hr_3(i) = HR(-38.7,3.84e4,3.65e1,-8.19e-2,1.14e-4,-6.96e-8,1.55e-11,Ts(i)); %J/mol
Hr_4(i) = HR(208.61,1.02e5,-9.61,-2.14e-1,3.2e-4,-1.84e-7,3.85e-11,Ts(i)); %J/mol
Hr_5(i) = HR(-206.1,3.17e4,1.69e2,-5.68e-1,6.42e-4,-3.3e-7,6.5e-11,Ts(i)); %J/mol
Hr_6(i) = HR(41.15,-3.21e4,2.15e2,-4.36e-1,4.37e-4,-2.16e-7,4.2e-11,Ts(i)); %J/mol
Hr_7(i) = HR(-247.31,6.38e4,-4.61e1,-1.32e-1,2.05e-4,-1.14e-7,2.3e-11,Ts(i)); %J/mol
Hr_8(i) = HR(-74.81,3.17e4,1.4e2,-4.83e-1,5.4e-4,-2.73e-7,5.28e-11,Ts(i)); %J/mol
Hr_9(i) = HR(131.3,-1e1,2.95e1,-8.49e-2,1.02e-4,-5.75e-8,1.23e-11,Ts(i)); %J/mol
Hr_10(i) = HR(173.3,3.21e4,-1.86e2,3.51e-1,-3.35e-4,1.59e-7,-2.97e-11,Ts(i)); %J/mol

end
     
%ode15s
%der=dvar/dt

%Fluid

der_a(17,1) = 0; %X
der_a(15,1) = 0; %CM_Ni
der_a(16,1) = 0; %CM_NiO

der_a(1,1) = -7.3e+08; %CH4
der_a(2,1) = 0; %H2
der_a(3,1) = 0; %CO
der_a(4,1) = 0; %CO2
der_a(5,1) = 0; %N2
der_a(6,1) = 0; %H2O
der_a(7,1) = 0; %Tf    
der_a(8,1)  = 0 ;%Ts
der_a(9,1)  = 0 ;%CH4
der_a(10,1) = 0 ;%H2
der_a(11,1) = 0 ;%CO
der_a(12,1) = 0 ;%CO2
der_a(13,1) = 0 ;%N2  
der_a(14,1) = 0 ;%H2O
for i = 2:n
%-------------------------------------------------------------------------------------------------------------------------------------    
%Fluid,by heterogeneous page 4-2. 2.1,equation 10 & 11 for line 378-390 the dimensional analysis check ok 2018/11/25

der_a(17,i) = (r1(i)+r2(i)+r3(i))/CM_NiO_in*Mw_Ni ;%X
der_a(15,i) = (r1(i)+r2(i)+r3(i))*Mw_Ni; %CM_Ni
der_a(16,i) = -(r1(i)+r2(i)+r3(i))*Mw_Ni; %CM_NiO
der_a(1,i) = ((-u_s(i)*(CM(1,i)-CM(1,i-1))/dz)-(k_mass(1,i)*a_specific*(CM(1,i-1)-CMs(1,i)))/porosity_bed) ;%CH4
der_a(2,i) = ((-u_s(i)*(CM(2,i)-CM(2,i-1))/dz)-(k_mass(2,i)*a_specific*CM_t(i)*(CM(2,i)-CMs(2,i)))/porosity_bed) ;%H2
der_a(3,i) = ((-u_s(i)*(CM(3,i)-CM(3,i-1))/dz)-(k_mass(3,i)*a_specific*CM_t(i)*(CM(3,i)-CMs(3,i)))/porosity_bed) ;%CO
der_a(4,i) = ((-u_s(i)*(CM(4,i)-CM(4,i-1))/dz)-(k_mass(4,i)*a_specific*CM_t(i)*(CM(4,i)-CMs(4,i)))/porosity_bed) ;%CO2
der_a(5,i) = 0 ;%N2
der_a(6,i) = ((-u_s(i)*(CM(6,i)-CM(6,i-1)))/dz-(k_mass(6,i)*a_specific*CM_t(i)*(CM(6,i)-CMs(6,i)))/porosity_bed) ;%H2O 
der_a(7,i) = ((-u_s(i)*(Tf(i)-Tf(i-1))/(CM_t(i)*Cp_t(i))/dz)+(hf(i)*a_specific*(Ts(i)-Tf(i))/(CM_t(i)*Cp_t(i))))/porosity_bed; %Tf    
%-------------------------------------------------------------------------------------------------------------------------------------
%Solid   
der_a(8,i) = ((V_cat*ro_solid*(r1(i)*(-Hr_1(i))+r2(i)*(-Hr_2(i))+r3(i)*(-Hr_3(i))+r4(i)*(-Hr_4(i))+r5(i)*(-Hr_5(i))+r6(i)*(-Hr_6(i))+r7(i)*(-Hr_7(i))+r8(i)*(-Hr_8(i))+r9(i)*(-Hr_9(i))))-(hf(i)*a_specific*(Ts(i)-Tf(i))))/(V_cat*Cps_t(i)*(1-porosity_catalyst)*ro_solid+V_cat*Cp_t(i)*porosity_catalyst*CM_t(i));%Ts
der_a(9,i) = (ro_solid*(-r3(i)-r4(i)-r6(i)-r7(i))) + k_mass(1,i)*a_specific*CM_t(i)*(CM(1,i-1)-CMs(1,i-1))/porosity_catalyst;%CH4
der_a(10,i) = (ro_solid*(-r1(i)+2*r3(i)+3*r4(i)+r5(i)+2*r6(i)+2*r7(i)+r8(i))) + k_mass(2,i)*a_specific*CM_t(i)*(CM(2,i-1)-CMs(2,i-1))/porosity_catalyst; %H2
der_a(11,i) = (ro_solid*(-r2(i)+r3(i)+r4(i)-r5(i)+2*r6(i)+r8(i)+2*r9(i))) + k_mass(3,i)*a_specific*CM_t(i)*(CM(3,i-1)-CMs(3,i-1))/porosity_catalyst; %CO
der_a(12,i) = (ro_solid*(-r2(i)+r5(i)-r6(i)-r9(i))) + k_mass(4,i)*a_specific*CM_t(i)*(CM(4,i-1)-CMs(4,i-1))/porosity_catalyst; %CO2
der_a(13,i) = 0; %N2
der_a(14,i) = (ro_solid*(r1(i)-r4(i)-r5(i)-r8(i))) + k_mass(6,i)*a_specific*CM_t(i)*(CM(6,i-1)-CMs(6,i-1))/porosity_catalyst; %H2O
end
der(1:17*n) = [der_a(1,1:n) der_a(2,1:n) der_a(3,1:n) der_a(4,1:n) der_a(5,1:n) der_a(6,1:n) der_a(7,1:n) der_a(8,1:n) der_a(9,1:n) der_a(10,1:n) der_a(11,1:n) der_a(12,1:n) der_a(13,1:n) der_a(14,1:n) der_a(15,1:n) der_a(16,1:n) der_a(17,1:n)];
der_c = der';
