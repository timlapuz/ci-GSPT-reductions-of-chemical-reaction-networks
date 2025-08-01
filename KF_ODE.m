% ODEs of the KF reaction
% 
% "Coordinate Independendent Model Reductions of Chemical
% Reaction Networks Based on Geometric Singular Perturbation Theory"
% T.E.F. Lapuz and M. Wechselberger 2025

function dydt = KF_ODE(t,y,alpha,beta,gamma,rho1,rho2,rho3,rho4,rho5,rho6)

% Renaming
x = y(1);
Y = y(2);
z = y(3);
s = y(4);

% The ODEs
N0 = [0; 0; inv(beta); 1];
f0 = [-beta*s*z + alpha*(1-s)];
F0 = N0*f0;
F1 = [rho1*s - rho2*x; rho3*x - rho4*Y; rho5*Y - rho6*z; gamma*(1-s)];
f = F0+F1;
dxdt = f(1);
dYdt = f(2);
dzdt = f(3);
dsdt = f(4);

% Concatenate into one vector
dydt = [dxdt; dYdt; dzdt; dsdt];
end