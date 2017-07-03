function coeffs = SOB(Y)
%%% Calculate stoich ratio for sulfur-oxidising bacteria
%%% using yield Y

syms HSm HCO3m H2O O2 SO4m2 CHO Hp

% H
r(1) = HSm + HCO3m + 2*H2O + 1.8*CHO + Hp == 0;
% O
r(2) = 3*HCO3m + 2*H2O + 2*O2 + 4*SO4m2 + 0.5*CHO  == 0;
% S
r(3) = HSm + SO4m2 == 0;
% C
r(4) = HCO3m + CHO == 0;
% e
r(5) = -HSm - HCO3m + Hp -2*SO4m2 == 0;

r(6) = HSm == -CHO * 24.6/Y;
r(7) = CHO == 1;

coeffs = struct2array(solve(r));

end