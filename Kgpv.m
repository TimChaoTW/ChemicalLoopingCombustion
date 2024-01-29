function Lambda = Kgpv(A,B,C,D,Tf)
% Gas phase thermal conductivity constants for species; form Paper Heterogenous page 6 table 5 check Ok but O2 unknow 2018/11/24
% Fluid vicosity, form  Paper Heterogenous page 5 table 4
% global TfK
Lambda = A*Tf^B / ( 1 + C/Tf + (D/Tf^2));