function Kcj = Kmtc(Sc)
% equilibrium constant & adsorption coefficient: model page10-table 5; supplementary page 17 table SI.5 & page 18 table SI.8
global ReK GK porosity_bed ro_t_massK
Kcj = (0.357*ReK^-0.359) * (Sc^-(2/3)) * (GK/(porosity_bed*ro_t_massK));