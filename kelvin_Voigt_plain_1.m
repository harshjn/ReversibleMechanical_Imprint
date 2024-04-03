clc
clear all
close all
% Parameters for the Kelvin-Voigt model
k = 50000;  % Elastic modulus in Pa
eta = 10000; % Viscosity in Pa.s

% Time parameters
dt = 1e-1; % Time step in seconds
t_mat = 0:dt:10; % Total time for simulation in seconds

% Total number of time steps
N = length(t_mat);

% Applied stress: Increasing linearly for the first half, then decreasing
halfN = round(N/2)-1;
%%

maxStrain=0.1
epsilon_mat = linspace(0,maxStrain,halfN+1); % Strain
%epsilon_mat(end+1) = epsilon_mat(halfN);

epsilon_dot_mat = diff(epsilon_mat(1:halfN+1))/dt; % Strain rate
%epsilon_dot_mat(end+1) =0;
epsilon_dot_old = -1*(epsilon_dot_mat(end));
epsilon_dot_mat(end+1) = -1*epsilon_dot_old;
%%
k_mat = k*epsilon_mat(1:halfN).^2;          % Corresponding values of k
eta_mat=eta*epsilon_mat(1:halfN);


% Create an interpolation functions
k_interp = @(x) interp1(epsilon_mat(1:halfN), k_mat, x, 'linear', 'extrap');
eta_interp = @(x) interp1(epsilon_mat(1:halfN),eta_mat,x,'linear','extrap');

% Initializing strain and strain rate arrays
sigma_mat = zeros(1,halfN);
%% Forward Cycle
% Simulation loop for the Kelvin-Voigt model

sigma_mat(1:halfN)=k_mat.*epsilon_mat(1:halfN)...
        + eta_mat.*epsilon_dot_mat(1:halfN);

plot(sigma_mat(1:halfN));
hold on; plot(sigma_mat(1:halfN)./epsilon_mat(1:halfN))

%% Return Cycle
tMat2 = [];
t2= 0;
%%
for i = 1: 200
epsilon0 = epsilon_mat(end);

eta_current = eta_interp(epsilon0);
k_current = k_interp(epsilon0);
DecayT = eta_current./k_current;

epsilon_relax = epsilon0*exp(-dt/DecayT);

epsilon_dot_relax = epsilon_relax-epsilon_mat(end);

sigma_relax=   k_current.*epsilon_relax...
        + eta_current.*epsilon_dot_relax;
epsilon_mat(end+1) = epsilon_relax;
sigma_mat(end+1) = sigma_relax;
epsilon_dot_mat(end+1) = epsilon_dot_relax;
t2 = t2+dt;
tMat2(end+1) = t2;
end

%%
plot(epsilon_mat(2:end),sigma_mat)
