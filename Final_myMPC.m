% MPC Control for SBMHS System
clear all; clc; close all

global N_pred N_con We Wdu n_st n_op n_ip
global xhat_k ek_f phy gama Lp_inf C_mat rk_pred uk_minus_1 ek_filt

% Load system variables and parameters
% load SBMHS_MIMO_ID_OE.mat
load Model_Data.mat
N_samples = 500; % Number of samples in simulation run

C_mat_plant = C_mat_OE; % Note: C matrix for plant and model are different
phy_plant = phy_OE;
gama_plant = gama_OE;

%Initialization of Simulation Related parameters
samp_T = 50; % Sampling time equals to 4 secs

% Initial state vector for simulating plant dynamics
[n_op, n_st_p] = size(C_mat_plant);
n_ip = length(Us);

%Create arrays for saving state variable time trajectories (Plant Simulation)
Xk = zeros(n_st_p, N_samples); % Array for saving state variable trajectories
Yk = zeros(n_op, N_samples); % will later add vk(measurement disturbance)

Uk(:,1) = Us;
Yk_abs(:,1) = Ys;

% Disturbance in the plant
Ds = Ys; % Initializing disturbace as the steady-state input
Dk = Ds - Ys;

%% Model Parameters and design (ARMAX MIMO Model)
C_mat = C_mat_ARMAX;
gama = gama_ARMAX;
phy = phy_ARMAX;
Lp_inf = K_mat_ARMAX;

[n_st, n_ip] = size(gama); % No. of state and input variables
[n_op, n_st] = size(C_mat); % No. of output

% Create arrays for saving state estimator time trajectories
dxk = zeros(n_st, N_samples); % Array for savinhg state variabels trajectories
duk = zeros(n_ip, N_samples);
ek = zeros(n_op, N_samples);
ek_f = zeros(n_op, N_samples); % Filtered Innovation Sequences


ip_noise_sigma_1 = std(duk_id(1,:)); % Measurement noise standard deviations
ip_noise_sigma_2 = std(duk_id(2,:));

%% Set MPC Tuning Parameters
alpha = 1;  beta = 1;  gamma = 0.9; delta = 0.9;
% wy11 = alpha / (meas_sigma_1 .^ 2);
% wy22 = beta / (meas_sigma_2 .^ 2);
wu11 = gamma / (ip_noise_sigma_1 .^ 2);
wu22 = delta / (ip_noise_sigma_2 .^2);
% 
% Wy = [wy11 0; 0 wy22];
Wdu = 10*[wu11 0; 0 wu22];  % Input weighting matrix
We = 30 * eye(n_op);        % Error weighting matrix
N_pred = 40;                % Prediction horizon
N_con = 5;                  % Control horizon
phy_e = 0.9*eye(n_op);      % Robustness Innovation filter parameter
uk_minus_1 = zeros(n_ip,1);

% Input constraints specifications
Ufk_min = zeros(N_con*n_ip,1);  % Future input vector
U_min = [0 0]';                 % Future input lower bound vector
u_min= U_min - Us;              % Future input move lower bound vector
U_max = [30 30]';               % Future input upper bound vector
u_max = U_max - Us;             % Future input move upper bound vector
delU = [2 2]';

Umin = [];
Umax = [];
delU_min = [];
delU_max = [];

for i = 1:N_con                 % Initialize input constraints vectors
    Umin = [Umin; u_min];
    delU_min = [delU_min; -delU];
    Umax = [Umax; u_max];
    delU_max = [delU_max; delU];
end

%% Setup matrices needed for including input move constraints
H_mat = zeros(n_ip, n_ip*N_con);
H_mat(1:n_ip, 1:n_ip) = eye(n_ip);
mat1 = diag(ones(n_ip*N_con,1));
mat2 = diag(ones(n_ip*(N_con-1),1), -n_ip);
delU_mat = mat1-mat2;
delU_mat0 = H_mat';
A_mat = [-delU_mat; delU_mat];

oldopts = optimset;
tolerance = 1e-8;
% options = optimset(oldopts,'MaxIter', 1e6, 'Display', 'iter', 'TolFun', tolerance, 'LargeScale', 'off', 'TolX', tolerance,'Algorithm','active-set');
options = optimset(oldopts,'MaxIter', 1e6, 'TolFun', tolerance, 'LargeScale', 'off', 'TolX', tolerance,'Algorithm','active-set');

% Setpoint trajectory related variables
rk = [0 0]';
rk_f = rk; % Filtered deviation setpoint trajectory
Rk(:,1) = rk_f + Ys;
phy_r = 0.9 * eye(n_op); % Setpoint filter parameter

%% Closed Loop Dynamics simulation of digital control system
for k = 2: N_samples
    k
    % MPC calculations at instant kT
    if k <= 50
        Dk = Ds - Ys;
        rk = [0 0]';
    elseif k > 50 && k < 300
        Dk = 0.9*Ds - Ys;
        rk = [0 0]';
    elseif k > 301
        Dk = 1.1*Ds - Ys;
        rk = [0 0]';
    end
    
    rk_f(:,k)=phy_r*rk_f(:,k-1) + (eye(n_op)-phy_r)*rk;
    Rk(:,k)=rk_f(:,k) + Ys; 
    vk = 0.5*randn(n_op,1);
 
        % Plant variables update
    Xk(:,k) = phy_plant * Xk(:,k-1) + gama_plant * duk(:,k-1); % samp_T has been used here
    Yk(:,k) = C_mat_plant * Xk(:,k) + Dk + vk; % This value gets looped for every iterations
    Yk_abs(:,k) = Yk(:,k) + Ys; % This abosolute value is only for plotting
    
    % Kalman predictor calculations for MPC over interval [kT, (k+1)T]
    dxk(:,k) = phy*dxk(:,k-1) + gama*duk(:,k-1) + Lp_inf*ek(:,k-1);
    
    ek(:,k) = Yk(:,k) - C_mat*dxk(:,k); % compute innovation 
    ek_f(:,k+1) = phy_e * ek_f(:,k) + (eye(n_op) - phy_e) * ek(:,k);
    
    ek_filt = ek_f(:,k);
    xhat_k = dxk(:,k);
    rk_pred = rk_f(:,k);
    
    if(k>2)
        uk_minus_1 = duk(:,k-1); 
    end
     
    vec1 = -delU_mat0 * uk_minus_1 - delU_min;
    vec2 = delU_mat0 * uk_minus_1 + delU_max;
    B_vec = [vec1; vec2];
    
    n1 = n_ip + 1; 
    n2 = (N_con - 1)*n_ip + 1;
    Ufk_0 = [Ufk_min(n1:end); Ufk_min(n2:end)];
    Ufk_min = fmincon('myMPC_ObjFn', Ufk_0, A_mat, B_vec, [], [], Umin, Umax, [], options);
    duk(:,k) = Ufk_min(1:n_ip);
    Uk(:,k) = Us + duk(:,k);

end   

figure, subplot(211),plot(Yk_abs(1,:),'k'),hold on
plot(Rk(1,:),'r'),hold on
xlabel('Sampling Instant(sec)'),ylabel('Temp_1(Kelvin)')
title('Plot of Plant Output Yk_1 Vs Setpoint_1');
legend('Yk_1', 'Setpoint Rk_1')
subplot(212), plot(Yk_abs(2,:),'k'), hold on
plot(Rk(2,:),'r'),hold on
xlabel('Sampling Instant(sec)'),ylabel('Temp_2(Kelvin)')
title('Plot of Plant Output Yk_2 Vs Setpoint_2');
legend('Yk_2', 'Setpoint Rk_2')


figure, subplot(211),stairs(Uk(1,(2:end)),'k'),grid
xlabel('Sampling Instant(sec)'),ylabel('Uk_1(Volts)')
title('Plot of Input Uk_1 and Uk_2');
subplot(212),stairs(Uk(2,(2:end)),'r'),grid
xlabel('Sampling Instant(sec)'),ylabel('Uk_2(Volts)')

figure, subplot(211),plot(ek(1),'k'),grid
xlabel('Sampling Instant(sec)'),ylabel('Uk_1(Volts)')
title('Plot of Input Uk_1 and Uk_2');
subplot(212),stairs(ek(2),'r'),grid
xlabel('Sampling Instant(sec)'),ylabel('Uk_2(Volts)')
