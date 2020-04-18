% Linear Quadratic Optimal Control (LQOC) program for SBMHS system
clear all; clc; close all

%% Plant Simulations
% Load SBMHS ARMAX and OE Model parameters and steady states points
load Model_Data.mat

N_samples = 500; % Number of samples in simulation run

phy_plant = phy_OE;
gama_plant = gama_OE;
C_mat_plant = C_mat_OE;


% Note: C matrix for plant and model are different
%Initialization of Simulation Related parameters
samp_T = 3; % Sampling time equals to 4 secs

% Initial state vector for simulating plant dynamics
[n_op_p, n_st_p] = size(C_mat_plant);
n_ip_p = length(Us);

%Create arrays for saving state variable time trajectories (Plant Simulation)
Xk_p = zeros(n_st_p, 1); % Array for saving state variable trajectories

Yk_p = C_mat_plant * Xk_p; % will later add vk(measurement disturbance)

% Disturbance in the plant
Ds = Ys; % Initializing disturbace as the steady-state input
Dk = Ds - Ys;

%% Model Parameters and design (ARMAX MIMO Model)
C_mat_model = C_mat_ARMAX;
gama_model = gama_ARMAX;
phy_model = phy_ARMAX;
Lp_inf = K_mat_ARMAX;

% Create dummy arrays for storing dynammic simulation results k'th column
% of these arrays corresponds to vector at k'th sampling instant
[n_st_m, n_ip_m] = size(gama_model); % No. of state and input variables
[n_op_m, n_st_m] = size(C_mat_model); % No. of output

% Create arrays for saving state estimator time trajectories
xk_hat_m = zeros(n_st_m, 1); % Array for savinhg state variabels trajectories

ek = Yk_p - C_mat_model * xk_hat_m;

%% Innovation Bias LQ Controller Design
% LQOC (Innovation Bias Formulation) Initialization
ek_f = zeros(n_op_m, 1); % Filtered Innovation Sequences

meas_sigma_1 = std(dyk_id(1,:));  % State noise standard deviations
meas_sigma_2 = std(dyk_id(2,:));

ip_noise_sigma_1 = std(duk_id(1,:)); % Measurement noise standard deviations
ip_noise_sigma_2 = std(duk_id(2,:));

%All these parameters > 0 are to choose to define relative importance of two outputs
alpha = 1;  beta = 1;  gama = 0.9; delta = 0.9;
wy11 = alpha / meas_sigma_1;
wy22 = beta / meas_sigma_2;
wu11 = gama / ip_noise_sigma_1;
wu22 = delta / ip_noise_sigma_2;

Wy = [wy11 0; 0 wy22];
Wu = 30*[wu11 0; 0 wu22]; % Input weighting matrix
Wx = 1*C_mat_model' * Wy * C_mat_model; % State weighting matrix

% Compute controller gain matrix by solving steady state Riccati equations
[G_inf, S_inf, EigVal_CL] = dlqr(phy_model, gama_model, Wx, Wu); % Compute controller gain matrix

% Robustness filter tuning parameters for LQOC
alpha_e = 0.95; 
alpha_r = 0.9; 
phy_e = alpha_e * eye(n_op_m); % Robustness Innovation filter parameter
phy_r = alpha_r * eye(n_op_m);

% Matrices required in Target state computation
Ku_mat = C_mat_model * inv(eye(n_st_m)-phy_model) * gama_model;
Ke_mat = C_mat_model * inv(eye(n_st_m)-phy_model) * Lp_inf + eye(n_op_m);

% % Computing target input and target states
% usk(:,1) = pinv(Ku_mat) * (rk_f(:,1) - Ke_mat * ek_f(:,1));
% xsk(:,1) = inv(eye(n_st_m) - phy_model)*(gama_model*usk(:,1) + Lp_inf*ek_f(:,1));

% uk_c(:,1) = usk(:,1) - G_inf * (xk_hat_m(:,1) - xsk(:,1));
uk_c = [0 0]';
Uk_c(:,1) = Us + uk_c;

U_min = [0 0]';
U_max = [30 30]';
u_min = U_min - Us;
u_max = U_max - Us;

% setpoint trajectory (Output is deviated in this case)
rk = [0 0]';
rk_f = rk;
Rk(:,1) = rk_f + Ys;

%% Closed loop system for Servo and Regulatory Simulations of the plant
for k = 1:N_samples
     if (k <= 50)
         Dk = (Ds -Ys);
         %Dk(:,k) = Ds - Ds;
         rk = [0 0]';
     elseif (k>50 && k<=270)
         Dk = 0.9*Ds -Ys;
         %Dk(:,k) = Ds - 0.9 * Ds;
         rk = [0 0]';
     elseif (k>270 && k<=N_samples)
         Dk = Ds -Ys;
         %Dk(:,k) = 1.1 * Ds - Ds;
         rk = [0 0]';
     end  
     
    % Filtered innovation for innovation bias implementation
    rk_f = phy_r * rk_f + (eye(n_op_p) - phy_r) * rk;
    
    Rk(:,k) = rk_f + Ys;
    
    vk = 0.5*randn(n_op_m,1);
    % Plant variables update
    Xk_p = phy_plant * Xk_p + gama_plant * uk_c;
    Yk_p = C_mat_plant * Xk_p + Dk + vk; % This value gets looped for every iterations
    Yk_p_abs(:,k) = Yk_p + Ys; % This abosolute value is only for plotting
    
    ek = Yk_p - C_mat_model * xk_hat_m; % Compute innovation
    ek_f = phy_e * ek_f + (eye(n_op_m) - phy_e) * ek;
    
    % Model variables update
    xk_hat_m = phy_model * xk_hat_m + gama_model * uk_c + Lp_inf * ek;
    
    
    % Control law calculations at instant kT
    % Computing target input and target states
    usk = pinv(Ku_mat) * (rk_f - Ke_mat * ek_f);
    xsk = inv(eye(n_st_m) - phy_model)*(gama_model*usk + Lp_inf*ek_f);
    uk_c = usk - G_inf * (xk_hat_m - xsk);
    
    % Impose Input bounds (constraints)
    if(uk_c(1)<= u_min(1)) % In SBMHS, minimum input voltage = 0 V
        uk_c(1) = u_min(1);
    elseif (uk_c(1) >= u_max(1)) % In SBMHS, maximum input voltage = 30 V
        uk_c(1) = u_max(1);
    end
    
    if(uk_c(2)<= u_min(2))
        uk_c(2) = u_min(2);
    elseif ( uk_c(2) >= u_max(2) )
        uk_c(2) = u_max(2);
    end
    
      % Input constraints handling
    Uk_c(:,k) = Us + uk_c; % Absolute vaule of input from the controller
    
end 

% Display simulations results gradually
figure, subplot(211), stairs(Uk_c(1,:)')
xlabel('Sampling Instant'), ylabel('U1(k)')
title('Plotting LQOC Inputs Uk_1 and Uk_2')
subplot(212), stairs(Uk_c(2,:)')
xlabel('Sampling Instant'), ylabel('U2(k)')

figure, plot(Rk(1,:), '-k'), hold on                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
plot(Yk_p_abs(1,:), '-r'), hold off
xlabel('Sampling Instant'), ylabel('Plot of R1(k) & Y1(k)')
title('Plottings of OE Plant output Yk_1 & Setpoint Rk_1')
legend('Setpoint Rk_1', 'Plant Output Yk_1')

figure, plot(Rk(2,:), '-k'), hold on
plot(Yk_p_abs(2,:), '-r'), hold off
xlabel('Sampling Instant'), ylabel('Plot of R2(k) & Y2(k)')
title('Plottings of OE Plant output Yk_2 & Setpoint Rk_2')
legend('Setpoint Rk_2', 'Plant Output Yk_2')