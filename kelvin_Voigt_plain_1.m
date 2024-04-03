clc
clear all
close all
% Parameters for the Kelvin-Voigt model
k = 50000;  % Elastic modulus in Pa
eta = 10000; % Viscosity in Pa.s

% Time parameters
dt = 1e-3; % Time step in seconds
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
n =(length(epsilon_mat)+1); 
t =t_mat(n); 
epsilon_current = epsilon_mat(end);
epsilon_dot_current= epsilon_dot_mat(end);
k_current =     k_interp(epsilon_current);
eta_current = eta_interp(epsilon_current);
%eta_dt = eta_current *
DecayT = eta_current./k_current;

epsilon_new = epsilon_current*exp(-dt/DecayT);
epsilon_dot_new = (epsilon_new -epsilon_current)/dt;

    while (epsilon_dot_new < epsilon_dot_old )
        epsilon_mat(end+1) = epsilon_mat(end) + dt*epsilon_dot_old;
        sigma_mat(end+1) = k_current*epsilon_mat(end)+ ...
            eta_current*epsilon_dot_old;
        epsilon_dot_mat(end+1) = epsilon_dot_old;

        n =(length(epsilon_mat)+1); 

        t =t_mat(n); 
        epsilon_current = epsilon_mat(end);
        epsilon_dot_current= epsilon_dot_mat(end);
        k_current =     k_interp(epsilon_current);
        eta_current = eta_interp(epsilon_current);
        %eta_dt = eta_current *
        DecayT = eta_current./k_current;

        epsilon_new = epsilon_current*exp(-dt/DecayT);
        epsilon_dot_new = (epsilon_new -epsilon_current)/dt;

end
%% Return cycle second half

n =(length(epsilon_mat)+1); 
t =t_mat(n); 
epsilon_current = epsilon_mat(end);
epsilon_dot_current= epsilon_dot_mat(end);
k_current =     k_interp(epsilon_current);
eta_current = eta_interp(epsilon_current);
%eta_dt = eta_current *
DecayT = eta_current./k_current;

epsilon_new = 0.1*exp(-dt/DecayT);
epsilon_dot_new = (epsilon_new -epsilon_current)/dt;





%% Plotting
close all
em=transpose(epsilon_mat(1:end-1));
sm=smooth(sigma_mat,1);
plot(em,sm)

title('Forward and return cycle')
xlabel('x')
ylabel('F')



%% plotting second half return cycle
%plot(em(halfN:end))


plot(em(1:halfN),smooth(sm(1:halfN),1000))
hold on
plot(em(halfN:end),smooth(sm(halfN:end),10))


%%
sm2=sm(5001:end)
em2=em(5001:end);
sm3(1:10)=linspace(0.699,0.3,10)
em3=[0.1*ones(10,1); em2];
sm4 = [sm3;  zeros(100,1)];
em4 = [em3;  linspace(em3(end),0,100)'];
