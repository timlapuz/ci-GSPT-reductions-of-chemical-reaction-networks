% ODEs of the irreversible MM.
% 
% "Coordinate Independendent Model Reductions of Chemical
% Reaction Networks Based on Geometric Singular Perturbation Theory"
% T.E.F. Lapuz and M. Wechselberger 2026

function dydt = MM_solver(t,y,alpha,beta,gamma)

% Renaming
s = y(1);
c = y(2);

% The ODEs
dsdt = (-s+(s+alpha)*c)*beta;
dcdt = s-(s+alpha+gamma)*c;

% Concatenate into one vector
dydt = [dsdt; dcdt];
end