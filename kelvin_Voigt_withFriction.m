clc
clear all
close all
% Parameters for the Kelvin-Voigt model
k = 50000;   % Elastic modulus in Pa
eta = 20000; % Viscosity in Pa.s

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
fric_mat = 4000.*epsilon_mat(1:halfN).^2;


% Create an interpolation functions
k_interp = @(x) interp1(epsilon_mat(1:halfN), k_mat, x, 'linear', 'extrap');
eta_interp = @(x) interp1(epsilon_mat(1:halfN),eta_mat,x,'linear','extrap');
fric_interp = @(x) interp1(epsilon_mat(1:halfN), fric_mat,x,'linear','extrap');

% Initializing strain and strain rate arrays
sigma_mat = zeros(1,halfN);
%% Forward Cycle
% Simulation loop for the Kelvin-Voigt model
close all
sigma_mat(1:halfN)=k_mat.*epsilon_mat(1:halfN)...
        + eta_mat.*epsilon_dot_mat(1:halfN)...
        -fric_mat.*epsilon_mat(1:halfN);

plot(epsilon_mat(1:halfN),sigma_mat(1:halfN));
hold on; plot(epsilon_mat(1:halfN),sigma_mat(1:halfN)./epsilon_mat(1:halfN))

%% Return Cycle
n =(length(epsilon_mat)+1); 
t =t_mat(n); 
epsilon_current = epsilon_mat(end);
epsilon_dot_current= epsilon_dot_mat(end);
k_current =     k_interp(epsilon_current);
eta_current = eta_interp(epsilon_current);
fric_current = fric_interp(epsilon_current);
%eta_dt = eta_current *
DecayT = eta_current./k_current;

epsilon_new = epsilon_current*exp(-dt/DecayT);
epsilon_dot_new = (epsilon_new -epsilon_current)/dt;

    while (epsilon_dot_new < epsilon_dot_old )
        epsilon_mat(end+1) = epsilon_mat(end) + dt*epsilon_dot_old;
        sigma_mat(end+1) = k_current*epsilon_mat(end)+ ...
            eta_current*epsilon_dot_mat(end)...
            - fric_current*epsilon_mat(end);
        epsilon_dot_mat(end+1) = epsilon_dot_old;

        n =(length(epsilon_mat)+1); 

        t =t_mat(n); 
        epsilon_current = epsilon_mat(end);
        epsilon_dot_current= epsilon_dot_mat(end);
        k_current =     k_interp(epsilon_current);
        eta_current = eta_interp(epsilon_current);
        %eta_dt = eta_current *
        DecayT = eta_current./k_current;

        epsilon_new = epsilon_current*exp(-dt/DecayT)
        epsilon_dot_new = (epsilon_new -epsilon_current)/dt

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

plot(epsilon_mat(5002:end),sigma_mat(5001:end))


%% Plotting
close all
em=transpose(epsilon_mat(1:end-1));
sm=sigma_mat;
%plot(em,sm./em)
%%
sm1 = sm(1:5000);
em1 = em(1:5000);

sm2=sm(5001:end);
em2=em(5001:end);
%sm3(1:10)=linspace(0.699,0.3,10)
%em3=[0.1*ones(10,1); em2];

tMat2 =linspace(0,5*DecayT,1000); 
em4=em2(end)*exp(-tMat2/DecayT);

dem4=(diff(em4)./diff(tMat2));
em4_ = em4(2:end);
sm4 = k_interp(em4_).*em4_+ ...
           eta_interp(em4_).*dem4...
           - fric_interp(em4_).*em4_;
% sm4 = -fric_interp(em4).*em4+k_interp(em4).*em4;
%+eta_interp(em4).*epsilon_dot_new;
em5=[em;em4'];
sm5=[sm sm4];
%%
plot(sm5)