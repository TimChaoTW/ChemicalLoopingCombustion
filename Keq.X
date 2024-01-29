function Kj = Keq(Kj0,DeltaH)
% equilibrium constant & adsorption coefficient: model page10-table 5; supplementary page 17 table SI.5 & page 18 table SI.8
% global TsK R_bar
global Rbar_Ts
% R_bar = 8.3151398; %Gas constant (Pa*m^3/(mol*K))
% Kj = Kj0*exp(-1*DeltaH/(R_bar*TsK));
Kj = Kj0*exp(-1*DeltaH/Rbar_Ts);