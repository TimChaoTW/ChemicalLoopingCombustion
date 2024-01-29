function D_AB = DC(M_AB,S_A,S_B)
% DIFFUSION COEFFICIENTS(m^2/s)
global Tdiff Pdiff
D_AB = Tdiff /(Pdiff*M_AB^0.5*(S_A^(1/3)+S_B^(1/3))^2); %DIFFUSION COEFFICIENTS(m^2/s)
