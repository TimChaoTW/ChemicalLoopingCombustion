function objfcn_Dynamic = Dynamic_PBR(dt,TotalTime)

global nt dti ro_bed porosity_bed d_p porosity_catalyst L n dz xyz
global CM_in yin_ch4 yin_h2 yin_co yin_co2 yin_N2 yin_h2o dtt Cps_t u_s  
global Pt_read  P_t time CM_CH4_in CM_H2_in CM_CO_in CM_CO2_in CM_N2_in CM_H2O_in Tf_in Pt_in_2 Tw_in_2 Tw_read T_w Ts_in 
global CMs_CH4_in CMs_H2_in CMs_CO_in CMs_CO2_in CMs_N2_in CMs_H2O_in CM_Ni_in CM_NiO_in

% L:Length of reactor (m); dz:Segment length
n = L/dz; %Number of segments

%Set inlet conditions
var0(1:n) = CM_CH4_in; %CH4 Fluid mole concentration
var0(n+1:2*n) = CM_H2_in; %H2 Fluid mole concentration
var0(2*n+1:3*n) = CM_CO_in; %CO Fluid mole concentration
var0(3*n+1:4*n) = CM_CO2_in; %CO2 Fluid mole concentration
var0(4*n+1:5*n) = CM_N2_in; %N2 Fluid mole concentration
var0(5*n+1:6*n) = CM_H2O_in; %H2O Fluid mole concentration
var0(6*n+1:7*n) = Tf_in; %Tf (K)
var0(7*n+1:8*n) =Ts_in; %Ts (K)
var0(8*n+1:9*n) = CMs_CH4_in; %CH4 Solid mole concentration
var0(9*n+1:10*n) = CMs_H2_in; %H2 Solid mole concentration
var0(10*n+1:11*n) =CMs_CO_in; %CO Solid mole concentration
var0(11*n+1:12*n) = CMs_CO2_in; %CO2 Solid mole concentration
var0(12*n+1:13*n) = CMs_N2_in; %N2 Solid mole concentration
var0(13*n+1:14*n) = CMs_H2O_in ; %H2O Solid mole concentration
var0(14*n+1:15*n) = CM_Ni_in; %Ni concentration
var0(15*n+1:16*n) = CM_NiO_in; %NiO concentration
var0(16*n+1:17*n) = 1e-6; %conversion of oxygen carrier

%Call ode15s for integration
odeopt=odeset('NonNegative', [1:17*n],'AbsTol',1e-3,'Stats','on');

tspan =[0 TotalTime]; % s, time of reactor:
time = TotalTime;
dtt =dt;

xyz = 0; % Debug counter --> view parameter ode15s' 't'

tic
[t, var] = ode15s (@SteadyState_PBR, tspan, var0,odeopt);
toc

%Profiles from integration
CM_CH4 = var(:,1:n);        %CH4 Fluid mole concentration mol/m^3
CM_H2 = var(:,n+1:2*n);     %H2 Fluid mole concentration mol/m^3
CM_CO = var(:,2*n+1:3*n);   %CO Fluid mole concentration mol/m^3
CM_CO2 = var(:,3*n+1:4*n);  %CO2 Fluid mole concentration mol/m^3
CM_N2 = var(:,4*n+1:5*n);   %N2 Fluid mole concentration mol/m^3
CM_H2O = var(:,5*n+1:6*n);  %H2O Fluid mole concentration mol/m^3
Tf = var(:,6*n+1:7*n); %Tf (K)
Ts = var(:,7*n+1:8*n); %Ts (K)
Cs_CH4 = var(:,8*n+1:9*n); %CH4 Solid mole fraction
Cs_H2 = var(:,9*n+1:10*n) ;%H2 Solid mole fraction
Cs_CO = var(:,10*n+1:11*n) ;%CO Solid mole fraction
Cs_CO2 = var(:,11*n+1:12*n); %CO2 Solid mole fraction
Cs_N2 = var(:,12*n+1:13*n) ;%O2 Solid mole fraction
Cs_H2O = var(:,13*n+1:14*n); %H2O Solid mole fraction
CM_Ni = var(:,14*n+1:15*n); %Ni concentration
CM_NiO = var(:,15*n+1:16*n); %NiO concentration
X = var(:,16*n+1:17*n); %conversion of oxygen carrier

% %Excel
% status=xlswrite('FluidMoleconcentration',CM_CH4,'CH4(mol)DividedBy(t)');
% status=xlswrite('FluidMoleconcentration',CM_H2,'H2(mol)DividedBy(t)');
% status=xlswrite('FluidMoleconcentration',CM_CO,'CO(mol)DividedBy(t)');
% status=xlswrite('FluidMoleconcentration',CM_CO2,'CO2(mol)DividedBy(t)');
% status=xlswrite('FluidMoleconcentration',CM_N2,'N2(mol)DividedBy(t)');
% status=xlswrite('FluidMoleconcentration',CM_H2O,'H2O(mol)DividedBy(t)');
% status=xlswrite('TemperaturePressure',Tf,'Tf(K)');
% status=xlswrite('TemperaturePressure',Pt,'Pt(Pa)');
% status=xlswrite('TemperaturePressure',Tw,'Tw(K)');
% status=xlswrite('TemperaturePressure',Ts,'Ts(K)');
% status=xlswrite('SolidMoleFraction',ys_CH4,'CH4');
% status=xlswrite('SolidMoleFraction',ys_H2,'H2');
% status=xlswrite('SolidMoleFraction',ys_CO,'CO');
% status=xlswrite('SolidMoleFraction',ys_CO2,'CO2');
% status=xlswrite('SolidMoleFraction',ys_N2,'N2');
% status=xlswrite('SolidMoleFraction',ys_H2O,'H2O');
% status=xlswrite('Niconcentration',CM_Ni,'CM_Ni');
% status=xlswrite('oxygencarrierconcentration',CM_NiO,'CM_NiO');
% status=xlswrite('conversion of oxygen carrier',X,'X');

%plot CH4 profiles
plot(t, CM_CH4(:,4),'g',t, CM_CH4(:,8),'b',t, CM_CH4(:,12),'c',t, CM_CH4(:,16),'m',t, CM_CH4(:,20),'y',t, CM_CH4(:,22),'k')
% title ('CH4 profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'1','png')
plot_sub('CH4 profiles','mol/m^3','CM_CH4')

%plot CM_H2 profiles
plot(t, CM_H2(:,4),'g',t, CM_H2(:,8),'b',t, CM_H2(:,12),'c',t, CM_H2(:,16),'m',t, CM_H2(:,20),'y',t, CM_H2(:,22),'k')
% title ('H2 profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'2','png')
plot_sub('H2 profiles','mol/m^3','CM_H2')

%plot CM_CO profiles
plot(t, CM_CO(:,4),'g',t, CM_CO(:,8),'b',t, CM_CO(:,12),'c',t, CM_CO(:,16),'m',t, CM_CO(:,20),'y',t, CM_CO(:,22),'k')
% title ('CO profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'3','png')
plot_sub('CO profiles','mol/m^3','CM_CO')

%plot CM_CO2 profiles
plot(t, CM_CO2(:,4),'g',t, CM_CO2(:,8),'b',t, CM_CO2(:,12),'c',t, CM_CO2(:,16),'m',t, CM_CO2(:,20),'y',t, CM_CO2(:,22),'k')
% title ('CO2 profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'4','png')
plot_sub('CO2 profiles','mol/m^3','CM_CO2')

%plot CM_N2 profiles
plot(t, CM_N2(:,4),'g',t, CM_N2(:,8),'b',t, CM_N2(:,12),'c',t, CM_N2(:,16),'m',t, CM_N2(:,20),'y',t, CM_N2(:,22),'k')
% title ('N2 profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'5','png')
plot_sub('N2 profiles','mol/m^3','CM_N2')

%plot CM_H2O profiles
plot(t, CM_H2O(:,4),'g',t, CM_H2O(:,8),'b',t, CM_H2O(:,12),'c',t, CM_H2O(:,16),'m',t, CM_H2O(:,20),'y',t, CM_H2O(:,22),'k')
title ('H2O profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'6','png')
plot_sub('H2O profiles','mol/m^3','CM_H2O')

%plot Tf profiles
plot(t, Tf(:,4),'g',t, Tf(:,8),'b',t, Tf(:,12),'c',t, Tf(:,16),'m',t, Tf(:,20),'y',t, Tf(:,22),'k')
% title ('Tf profiles')
% xlabel ('Time sec')
% ylabel ('K')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'7','png')
plot_sub('Tf profiles','K','Tf')

%plot Ts profiles
plot(t, Ts(:,5),'g',t, Ts(:,8),'b',t, Ts(:,12),'c',t, Ts(:,16),'m',t, Ts(:,20),'y',t, Ts(:,22),'k')
% title ('Ts profiles')
% xlabel ('Time sec')
% ylabel ('K')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'8','png')
plot_sub('Ts profiles','K','Ts')

%plot ys_CH4 profiles
plot(t, Cs_CH4(:,4),'g',t, Cs_CH4(:,8),'b',t, Cs_CH4(:,12),'c',t, Cs_CH4(:,16),'m',t, Cs_CH4(:,20),'y',t, Cs_CH4(:,22),'k')
% title ('ys CH4 profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'9','png')
plot_sub('Cs CH4 profiles','','Cs_CH4')

%plot ys_H2 profiles
plot(t, Cs_H2(:,4),'g',t, Cs_H2(:,8),'b',t, Cs_H2(:,12),'c',t, Cs_H2(:,16),'m',t, Cs_H2(:,20),'y',t, Cs_H2(:,22),'k')
% title ('ys H2 profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'10','png')
plot_sub('Cs H2 profiles','','Cs_H2')

%plot ys_CO profiles
plot(t, Cs_CO(:,4),'g',t, Cs_CO(:,8),'b',t, Cs_CO(:,12),'c',t, Cs_CO(:,16),'m',t, Cs_CO(:,20),'y',t, Cs_CO(:,22),'k')
% title ('ys CO profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'11','png')
plot_sub('Cs CO profiles','','Cs_CO')

%plot ys_CO2 profiles
plot(t, Cs_CO2(:,4),'g',t, Cs_CO2(:,8),'b',t, Cs_CO2(:,12),'c',t, Cs_CO2(:,16),'m',t, Cs_CO2(:,20),'y',t, Cs_CO2(:,22),'k')
% title ('ys CO2 profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'12','png')
plot_sub('Cs CO2 profiles','','Cs_CO2')

%plot ys_N2 profiles
plot(t, Cs_N2(:,4),'g',t, Cs_N2(:,8),'b',t, Cs_N2(:,12),'c',t, Cs_N2(:,16),'m',t, Cs_N2(:,20),'y',t, Cs_N2(:,22),'k')
title ('Cs N2 profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'13','png')
plot_sub('Cs N2 profiles','','Cs_N2')

%plot ys_H2O profiles
plot(t, Cs_H2O(:,4),'g',t, Cs_H2O(:,8),'b',t, Cs_H2O(:,12),'c',t, Cs_H2O(:,16),'m',t, Cs_H2O(:,20),'y',t, Cs_H2O(:,22),'k')
% title ('ys H2O profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'14','png')
plot_sub('Cs H2O profiles','','Cs_H2O')

%plot CM_Ni profiles
plot(t, CM_Ni(:,4),'g',t, CM_Ni(:,8),'b',t, CM_Ni(:,12),'c',t, CM_Ni(:,16),'m',t, CM_Ni(:,20),'y',t, CM_Ni(:,22),'k')
% title ('CM Ni profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'15','png')
plot_sub('Ni profiles','mol/m^3','CM_Ni')

%plot CM_NiO profiles
plot(t, CM_NiO(:,4),'g',t, CM_NiO(:,8),'b',t, CM_NiO(:,12),'c',t, CM_NiO(:,16),'m',t, CM_NiO(:,20),'y',t, CM_NiO(:,22),'k')
% title ('CM NiO profiles')
% xlabel ('Time sec')
% ylabel ('mol/m^3')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'16','png')
plot_sub('NiO profiles','mol/m^3','CM_NiO')

%plot X profiles
plot(t, X(:,4),'g',t, X(:,8),'b',t, X(:,12),'c',t, X(:,16),'m',t, X(:,20),'y',t, X(:,22),'k')
% title ('X profiles')
% xlabel ('Time sec')
% legend('Z=2m','Z=4m','Z=6m','Z=8m','Z=10m','Z=11m')
% saveas(gcf,'17','png')
plot_sub('X profiles','','X')

% plot total ys
% plot(t, ys_CH4(:,22),'r',t,ys_H2(:,22),'b',t, ys_H2O(:,22),'c',t, ys_CO(:,22),'m',t, ys_CO2(:,22),'y')
% title ('gas fraction')
% xlabel ('Time sec')
% legend('CH4','H2','H2O','CO','CO2')
% saveas(gcf,'Mixer','png')

%Objective functions
objfcn_Dynamic(:,1) = CM_CO2(:,n);

%for i = 1:time/dt+1
%objfcn_Dynamic(i,2) = (FM_CH4(i,1)-FM_CH4(i,n+1))/FM_CH4(i,1);
%end
