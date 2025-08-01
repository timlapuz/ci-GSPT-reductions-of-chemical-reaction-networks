% Generates Figure 2b of "Coordinate Independendent Model Reductions of Chemical
% Reaction Networks Based on Geometric Singular Perturbation Theory"
% T.E.F. Lapuz and M. Wechselberger 2025

%% Setting up
% Parameter values
alpha = 0.75;
beta = 1;
gamma = 0.005;

% Final time
t_end = 3000;

% Integration time
tspan = [0 t_end];

%% Full 2D system
% ICs
y0 = [1,0];

% Integration
[t,y] = ode15s(@(t,y) MM_ODE(t,y, alpha, beta, gamma), tspan, y0);

% Plotting phase plane
figure(1); 
hold on;
plot(y(:,1),y(:,2),'Color',[0.2 0.2 0.8],'LineWidth',4); 
set(gca,'fontsize', 16) 
xlabel('$s$','Interpreter','Latex', 'FontSize', 20);  
ylabel('$c$','Interpreter','Latex', 'FontSize', 20);
grid on; 
xlim([0,1])
ylim([0,0.5])

%% ci-GSPT approximation
% Calculation of the IC on the critical manifold
delta = alpha + beta - 1;
s0_CM = (-delta + sqrt(delta^2+4*alpha))/2;
c0_CM = s0_CM/(s0_CM+alpha);

% Integration time
tspan = [0 t_end];

% ICs
y0 = [s0_CM];

% Integration
[t,y] = ode15s(@(t,y) MM_ODE_GSPT(t,y, alpha, beta, gamma), tspan, y0);

% Define
sp = y(:,1);
cp = (sp./(sp+alpha));

% Phase plane plots
figure(1)
plot(sp,cp,'Color',[0.2 0.8 0.2],'LineStyle','--','LineWidth',4); 
legend('Full', 'GSPT', 'FontSize', 14)

