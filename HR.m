function DeltaH = HR(Redu,a0,a1,a2,a3,a4,a5,Ts)
% equilibrium constant & adsorption coefficient: model page10-table 5; supplementary page 17 table SI.5 & page 18 table SI.8
% global TsK
RT = 298;
DeltaH = (a0*(Ts-RT)+a1/2*(Ts^2-RT^2)+a2/3*(Ts^3-RT^3)+a3/4*(Ts^4-RT^4)+a4/5*(Ts^5-RT^5)+a5/6*(Ts^6-RT^6))*1e-3 + Redu*1e3; %J/mol