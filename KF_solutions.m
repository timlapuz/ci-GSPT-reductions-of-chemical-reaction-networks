% Generates Figure 3b of "Coordinate Independendent Model Reductions of Chemical
% Reaction Networks Based on Geometric Singular Perturbation Theory"
% T.E.F. Lapuz and M. Wechselberger 2025

%% Setting up
% Parameter values
alpha =  4e-3;
beta = 1;
rho1 = 5e-6;
rho2 = 1e-6;
rho3 = 1e-5;
rho4 = 1e-6;
rho5 = 1e-6;
rho6 = 1e-6;
gamma = rho6;
gamma_perturb = 0.5*gamma;
gamma = gamma + gamma_perturb; % comment out for Figure 3a



% Integration time
t_end = 100000000;
tspan = [0 t_end];

% Specify error tolerance for integration step
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

z0 = (-(alpha+1-beta)+sqrt((alpha+1-beta)^2+4*alpha))/(2*beta);
s0 = alpha/(alpha+beta*z0);

 
% ICs
y0 = [1;1;z0;s0];

%% Full 4D system

% ICs
z0 = (-(alpha+1-beta)+sqrt((alpha+1-beta)^2+4*alpha))/(2*beta);
s0 = alpha/(alpha+beta*z0);
y0 = [1;1;z0;s0];

% Integration
[t,y] = ode15s(@(t,y) KF_ODE(t,y, alpha,beta, gamma,rho1,rho2,rho3,rho4,rho5,rho6), tspan, y0, opts);

% Rename
x = y(:,1);
Y = y(:,2);
z = y(:,3);
s = y(:,4);

% Plot figure
figure(1)
plot(t,z,'LineWidth',4,'Color',[0.2 0.2 0.8],'LineStyle','-')
hold on;
grid on;


%% ci-GSPT approximation
% ICs
y0 = y0(1:3);

[t,y] = ode15s(@(t,y) KF_ODE_GSPT(t,y, alpha,beta, gamma, rho1,rho2,rho3,rho4,rho5,rho6), tspan, y0, opts);

% Rename
x = y(:,1);
Y = y(:,2);
z = y(:,3);

% Plot figure
figure(1)
plot(t,z,'LineWidth',4,'Color',[0.2 0.8 0.2],'LineStyle','--'); 
hold on;
grid on;
legend('Full $z$','GSPT $z$','interpreter','latex')
set(gca,'FontSize',17)
xlim([4e7,6e7])

