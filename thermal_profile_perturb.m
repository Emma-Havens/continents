function perturb = thermal_profile_perturb(T)

H       =    350;           % kJ / mol
R       =    0.008314;      % kJ / (K * mol)
T_t     =    285;           % K (dimensional top TBL)
Del_T   =    3215;          % K (dimensional)
T_dim   =    T*Del_T+T_t;   % K (dimensional T)
C       =    14.73161684;

MU_0    =    10E20 * exp(-H/(C*R*T_t));
interm  =    exp(H/(C*R*T_dim));
MU_T    =    MU_0 * interm;

perturb =    MU_T / 10E18;  % dimensionless

end

% T_star  =    5611.314028*(T^2) - 1229.834638*(T) + 3421.854812;
% perturb =    exp((-H / (R * T_star)) + 10);


