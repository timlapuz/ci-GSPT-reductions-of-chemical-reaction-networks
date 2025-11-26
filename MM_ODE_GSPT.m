% ODE of a ci-GSPT reduction (gamma<<1,alpha,beta=O(1)) of the irreversible MM
% 
% "Coordinate Independendent Model Reductions of Chemical
% Reaction Networks Based on Geometric Singular Perturbation Theory"
% T.E.F. Lapuz and M. Wechselberger 2026


function dydt = MM_ODE_GSPT(t,y,alpha,beta,gamma)

% Renaming
s = y(1);

% The ODEs
dsdt = -(gamma*beta*s*(alpha + s))/(alpha*beta + (alpha + s)^2);

% Concatenate into one vector
dydt = [dsdt];
end