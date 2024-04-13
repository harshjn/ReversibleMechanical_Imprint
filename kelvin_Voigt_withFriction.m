clc
clear all
close all
% Parameters for the Kelvin-Voigt model
k = 2*50000;   % Elastic modulus in Pa
eta = 10000; % Viscosity in Pa.s

% Time parameters
dt = 1e-1; % Time step in seconds
t_mat = 0:dt:10; % Total time for simulation in seconds

% Total number of time steps
N = length(t_mat);

% Applied strain: Increasing linearly for the first half. Second half: Applied stress set to zero, exponential decay in strain
halfN = round(N/2)-1;
%%

maxStrain=0.1;
x_mat = linspace(0,maxStrain,halfN+1); % Strain
%x_mat(end+1) = x_mat(halfN);

x_dot_mat = diff(x_mat(1:halfN+1))/dt; % Strain rate
%x_dot_mat(end+1) =0;
x_dot_old = -1*(x_dot_mat(end));
x_dot_mat(end+1) = -1*x_dot_old;
%%
k_mat = k*x_mat(1:halfN).^2;          % Corresponding values of k
eta_mat=eta*x_mat(1:halfN);
fric_mat = 4000.*x_mat(1:halfN).^2;


% Create an interpolation functions
k_interp = @(x) interp1(x_mat(1:halfN), k_mat, x, 'linear', 'extrap');
eta_interp = @(x) interp1(x_mat(1:halfN),eta_mat,x,'linear','extrap');
fric_interp = @(x) interp1(x_mat(1:halfN), fric_mat,x,'linear','extrap');

% Initializing strain and strain rate arrays
f_mat = zeros(1,halfN);
%% Forward Cycle
% Simulation loop for the Kelvin-Voigt model
close all
f_mat(1:halfN)=k_mat.*x_mat(1:halfN)...
        + eta_mat.*x_dot_mat(1:halfN)...
        +fric_mat;

plot(x_mat(1:halfN),f_mat(1:halfN));
hold on; plot(x_mat(1:halfN),f_mat(1:halfN)./x_mat(1:halfN))

%% Return Cycle
tMat2 = [];
t2= 0;
%%
for i = 1: 200
x0 = x_mat(end);

eta_eff  = eta_interp(x0);
k_eff    = k_interp(x0);
fric_eff = fric_interp(x0);
DecayT = (eta_eff)./(k_eff-fric_eff/x0);

x_relax = x0*exp(-dt/DecayT);

x_dot_relax = x_relax-x_mat(end);

f_relax   =   k_eff.*x_relax...
                + eta_eff.*x_dot_relax...
                - fric_eff;
x_mat(end+1) = x_relax;
f_mat(end+1) = f_relax;
x_dot_mat(end+1) = x_dot_relax;
t2 = t2+dt;
tMat2(end+1) = t2;
end

plot(x_mat)
%%
plot(x_mat(2:end),f_mat,'o--')
