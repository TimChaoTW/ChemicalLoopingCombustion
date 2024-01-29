function Cpi = Kgphc(A,B,C,D,E,Tf)
% Gas phase heat capacity constants; From dynamic page 256 Table 12.3; check OK 2018/11/24
% global TfK
Cpi = A + B*(C/(Tf*sinh((C/Tf)^2))) + D*(E/(Tf*sinh((E/Tf)^2)))*1e-3; %J/(mol*K)
