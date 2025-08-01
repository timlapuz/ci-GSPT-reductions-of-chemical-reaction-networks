% ODE of a ci-GSPT reduction (gamma<<1,alpha,beta=O(1)) of the KF reaction.
% 
% "Coordinate Independendent Model Reductions of Chemical
% Reaction Networks Based on Geometric Singular Perturbation Theory"
% T.E.F. Lapuz and M. Wechselberger 2025

function dydt = KF_ODE_GSPT(t,y,alpha,beta,gamma,rho1,rho2,rho3,rho4,rho5,rho6)

% Renaming
x = y(1);
Y = y(2);
z = y(3);

% The ODEs
K = (alpha+beta*z)/((alpha+beta*z)^2+alpha);
dxdt = alpha*rho1/(alpha+beta*z)-rho2*x;
dYdt = rho3*x-rho4*Y;
dzdt = K*((alpha+beta*z)*(rho5*Y-rho6*z) - gamma*z);

% Concatenate into one vector
dydt = [dxdt; dYdt; dzdt];
end