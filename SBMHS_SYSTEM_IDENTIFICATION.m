%% Development of a MIM ARMAX Model and OE Model from perturbation data
clear all; clc; close all

% Load excel file data
sbmhs_data = xlsread('sbmhs.xlsx');
[row, col] = size(sbmhs_data);

% Graphical Represenation of measurement data and manipulated inputs
figure, subplot(221), plot(sbmhs_data(:,7)), 
xlabel('Sampling Instants'), ylabel('Temperature (deg. Celcius)')
title('SBMHS Measurement data: Output Y1(k)')
subplot(223), plot(sbmhs_data(:,3)),
xlabel('Sampling Instants'), ylabel('Amplitude = Us + 2')
title('SBMHS Manipulated Inputs data: Input u1(k)')
subplot(222), plot(sbmhs_data(:,9)),
xlabel('Sampling Instants'), ylabel('Temperature (deg. Celcius)')
title('SBMHS Measurement data: Output Y2(k)')
subplot(224), plot(sbmhs_data(:,5)),
xlabel('Sampling Instants'), ylabel('Amplitude = Us + 3')
title('SBMHS Manipulated inputs data: Input u2(k)')

%% System Identification
sbmhs_ID_data = sbmhs_data(1:(0.75*row), :); % 75% of data taken for system identification
[row_id, col_id] = size(sbmhs_ID_data);

N_samples_id = row_id;
samp_T = 3; % in seconds

% Measurement/ Output data
Yk1_id = sbmhs_ID_data(:,7);
Yk2_id = sbmhs_ID_data(:,9);
Yk_id = [Yk1_id(309:end)' ; Yk2_id(309:end)'];

% Input data
Uk1_id = sbmhs_ID_data(:,3);     
Uk2_id = sbmhs_ID_data(:,5);     
Uk_id = [Uk1_id(309:end)';Uk2_id(309:end)'];

% Steady state operating points
Ys1 = mean(Yk1_id(309-50:309));   % Taking average of end 50 samples of Yk1
Ys2 = mean(Yk2_id(309-50:309));   % Taking average of end 50 samples of Yk2
Ys = [Ys1; Ys2];                  % Steady state output operating point

Us1 = 10;
Us2 = 10;
Us = [Us1; Us2];                % Steady state input operating point

[n_op, N_samples_id] = size(Yk_id);
[n_ip, N_samples_id] = size(Uk_id);

%% Generate perturbation data for model identification
for k = 1:N_samples_id
   duk_id(:,k) = Uk_id(:,k) - Us;  %%% confusion in taking steady state points
   dyk_id(:,k) = Yk_id(:,k) - Ys;
end

%% Create a data object according SysId Toolbox format
% Data object for MISO-I
id_data_1 = iddata(dyk_id(1,:)', duk_id', samp_T);

% Data object for MISO-II
id_data_2 = iddata(dyk_id(2,:)', duk_id', samp_T);

% Data object for MIMO
id_data = iddata(dyk_id', duk_id', samp_T);


%% Identify Model order for MISO ARMAX models
% Polynomial order for ARMAX and OE Model for MISO-I
nb_1_armax = [4 4];     % Order of B polynomial(numerator) + 1 (Ny-by-Nu)
nb_1_oe = [2 2];
nf_1 = [2 2];
na_1 = 4;         % Oder of F polynomial(denominator) (Ny-by-Nu)
nc_1 = 2;         % Order of C-polynomial
nk_1 = [2 5];     % Time delays w.r.t each input (Ny-by-Nu)

% Polynomial order for ARMAX and OE Model for MISO-II
nb_2_armax = [4 4];     % Order of B polynomial(numerator) + 1 (Ny-by-Nu)
nb_2_oe = [2 2];
nf_2 = [2 2];
na_2 = 4;         % Oder of F polynomial(denominator) (Ny-by-Nu)
nc_2 = 2;         % Order of C-polynomial
nk_2 = [7 3];     % Time delays w.r.t each input (Ny-by-Nu)

% Armax Model for MISO-1
ARMAX_model_1 = armax(id_data_1, [na_1 nb_1_armax nc_1 nk_1]);

% Armax Model for MISO-II
ARMAX_model_2 = armax(id_data_2, [na_2 nb_2_armax nc_2 nk_2]);

% OE Model for MISO-1
OE_model_1 = oe(id_data_1, [nb_1_oe nf_1 nk_1]);

% OE Model for MISO-II
OE_model_2 = oe(id_data_2, [nb_2_oe nf_2 nk_2]);

% Predictions for MISO-1 on identification data
dyk1_hat_armax= predict(ARMAX_model_1, id_data_1, 1);

dyk1_hat_oe = predict(OE_model_1, id_data_1,1);


% Predictions for MISO-II on identification data
dyk2_hat_armax= predict(ARMAX_model_2, id_data_2, 1);

dyk2_hat_oe = predict(OE_model_2, id_data_2,1);


figure, subplot(211), plot(dyk_id(1,:), 'xr')
hold on, plot(dyk1_hat_armax.y,'-k')
plot(dyk1_hat_oe.y, '-b'), hold off
xlabel('Sampling Instants'), ylabel('Temperature (deg. Celcius)')
title('SBMHS System Identification: MISO-I')
legend('Measurement', 'ARMAX Prediction', 'OE Prediction')

subplot(212), plot(dyk_id(2,:), 'xr')
hold on, plot(dyk2_hat_armax.y,'-k')
plot(dyk2_hat_oe.y, '-b'), hold off
xlabel('Sampling Instants'), ylabel('Temperature (deg. Celcius)')
title('SBMHS System Identification: MISO-II')
legend('Measurement', 'ARMAX Prediction', 'OE Prediction')

% Residuals auto/cross correlation for MISO-1 ARMAX Model
figure, subplot(221), resid(ARMAX_model_1, id_data_1, 'CORR');
% Residuals auto/cross correlation for MISO-1 OE Model
subplot(222), resid(OE_model_1, id_data_1, 'CORR');

% Residuals auto/cross correlation for MISO-II ARMAX Model
subplot(223), resid(ARMAX_model_2, id_data_2, 'CORR');
% Residuals auto/cross correlation for MISO-II OE Model
subplot(224), resid(OE_model_2, id_data_2, 'CORR');

%% MISO-I - ARMAX & OE MODEL
% Convert ARMAX model to state space form for MISO-I
[phy1_ARMAX, gama1_ARMAX, C_mat1_ARMAX, D_mat1_ARMAX, K_mat1_ARMAX, X01_ARMAX] = ssdata(ARMAX_model_1);

% Convert OE model to state space form for MISO-I
[phy1_OE, gama1_OE, C_mat1_OE, D_mat1_OE, X01_OE] = ssdata(OE_model_1);

% Convert to Matlab state space and transfer function objects MISO-I/ ARMAX
dss1_ARMAX = ss(phy1_ARMAX, gama1_ARMAX, C_mat1_ARMAX, D_mat1_ARMAX, samp_T);
dtf1_ARMAX = tf(dss1_ARMAX);

% Convert to Matlab state space and transfer function objects MISO-I/ OE
dss1_OE = ss(phy1_OE, gama1_OE, C_mat1_OE, D_mat1_OE, samp_T);
dtf1_OE = tf(dss1_OE);

%% MISO-II - ARMAX AND OE MODEL
% Convert ARMAX model to state space form for MISO-II
[phy2_ARMAX, gama2_ARMAX, C_mat2_ARMAX, D_mat2_ARMAX, K_mat2_ARMAX, X02_ARMAX] = ssdata(ARMAX_model_2);

% Convert OE model to state space form for MISO-II
[phy2_OE, gama2_OE, C_mat2_OE, D_mat2_OE, X02_OE] = ssdata(OE_model_2);

% Convert to Matlab state space and transfer function objects MISO-II/ ARMAX
dss2_ARMAX = ss(phy2_ARMAX, gama2_ARMAX, C_mat2_ARMAX, D_mat2_ARMAX, samp_T);
dtf2_ARMAX = tf(dss2_ARMAX);

% Convert to Matlab state space and transfer function objects MISO-II/ OE
dss2_OE = ss(phy2_OE, gama2_OE, C_mat2_OE, D_mat2_OE, samp_T);
dtf1_OE = tf(dss2_OE);

%% MIMO state realization for ARMAX Model
zeros_mat = zeros(4);
zeros_vec = zeros(1,4);
phy_ARMAX = [phy1_ARMAX zeros_mat; zeros_mat phy2_ARMAX];
gama_ARMAX = [gama1_ARMAX; gama2_ARMAX];
C_mat_ARMAX = [C_mat1_ARMAX zeros_vec; zeros_vec C_mat2_ARMAX];
D_mat_ARMAX = [D_mat1_ARMAX; D_mat2_ARMAX];
K_mat_ARMAX = [K_mat1_ARMAX zeros_vec'; zeros_vec' K_mat2_ARMAX];

% Convert to Matlab state space and transfer function objects for MIMO
dss_ARMAX = ss(phy_ARMAX, gama_ARMAX, C_mat_ARMAX, D_mat_ARMAX, samp_T);
dtf_ARMAX = tf(dss_ARMAX);


%% MIMO state realization for OE Model
zeros_mat = zeros(4);           %%%%%%% Here why size 24*24??? Why not 12*12
zeros_vec = zeros(1,4);
phy_OE = [phy1_OE zeros_mat; zeros_mat phy2_OE];
gama_OE = [gama1_OE; gama2_OE];
C_mat_OE = [C_mat1_OE zeros_vec; zeros_vec C_mat2_OE];
D_mat_OE = [D_mat1_OE; D_mat2_OE];

% Convert to Matlab state space and transfer function objects for MIMO
dss_OE = ss(phy_OE, gama_OE, C_mat_OE, D_mat_OE, samp_T);
dtf_OE = tf(dss_OE);


%% Model Validation
sbmhs_Val_data = sbmhs_data((row_id+1):end, :); % Remaining 25% data taken for validation
[row_val, col_val] = size(sbmhs_Val_data);

N_samples_val = row_val;

% Measurement/ Output data
Yk1_val = sbmhs_Val_data(:,7);
Yk2_val = sbmhs_Val_data(:,9);
Yk_val = [Yk1_val' ; Yk2_val'];

% Input data
Uk1_val = sbmhs_Val_data(:,3);
Uk2_val = sbmhs_Val_data(:,5);
Uk_val = [Uk1_val'; Uk2_val'];

% Steady state operating points will be same as SysID data

[n_op, N_samples_val] = size(Yk_val);
[n_ip, N_samples_val] = size(Uk_val);

% Generate perturbation data for model validation
for k = 1:N_samples_val
    duk_val(:,k) = Uk_val(:,k) - Us;
    dyk_val(:,k) = Yk_val(:,k) - Ys;
end

%% MISO-I Predictions on Validation data for ARMAX Model and OE Model
% Infinite horizon predictions on validation data
n_st_ARMAX = length(phy1_ARMAX);
n_st_OE = length(phy1_OE);
xk_ARMAX = zeros(n_st_ARMAX, N_samples_val);
xk_OE = zeros(n_st_OE, N_samples_val);

% Predictions for validation for MISO-I ARMAX Model
yk1_hat_ARMAX(1) = C_mat1_ARMAX * xk_ARMAX(:,1);
for k = 1:N_samples_val-1
   xk_ARMAX(:,k+1) = phy1_ARMAX  * xk_ARMAX(:,k) + gama1_ARMAX * duk_val(:,k);
   yk1_hat_ARMAX(k+1) = C_mat1_ARMAX * xk_ARMAX(:,k+1);
end
yk1_hat_inf_ARMAX = yk1_hat_ARMAX;

% ARMAX Model: One Step Prediction on validation data for MISO-I ARMAX
xk_ARMAX = zeros(n_st_ARMAX, N_samples_val);
yk1_hat_ARMAX(1) = C_mat1_ARMAX * xk_ARMAX(:,1);
for k = 1: N_samples_val-1
    err_1(k) = dyk_val(1,k) - C_mat1_ARMAX * xk_ARMAX(:,k);
    xk_ARMAX(:,k+1) = phy1_ARMAX * xk_ARMAX(:,k) + gama1_ARMAX * duk_val(:,k) + K_mat1_ARMAX * err_1(k);
    yk1_hat_ARMAX(k+1) = C_mat1_ARMAX * xk_ARMAX(:,k+1) + err_1(k);
end
yk1_hat_one_ARMAX = yk1_hat_ARMAX;

% Predictions for validation for MISO-I OE Model
xk_OE = zeros(n_st_OE, N_samples_val);
yk1_hat_OE(1) = C_mat1_OE * xk_OE(:,1);
for k = 1: N_samples_val-1
    xk_OE(:,k+1) = phy1_OE * xk_OE(:,k) + gama1_OE * duk_val(:,k);
    yk1_hat_OE(k+1) = C_mat1_OE * xk_OE(:,k+1);
end

% Graphical visualization and comaparision for MISO-I
figure, subplot(2,2,[1 2]), plot(dyk_val(1,:),'xr'),grid
hold on, plot(yk1_hat_inf_ARMAX,'-k'), grid
plot(yk1_hat_one_ARMAX,'-b'), grid
plot(yk1_hat_OE, 'g'), hold off
legend('Measure Output','Infinite Horizon Prediction','One Step Prediction','OE Prediction');
xlabel('Sampling Instant'), ylabel('Measured and Predicted Outputs')
title('Model Validation: MISO-I Predictions Comparision')

subplot(223), stairs(duk_id(1,:)')
xlabel('Sampling Instant'),ylabel('u_1(k)')
title('Model Validation: Manipulated Inputs');
subplot(224), stairs(duk_id(2,:)')
xlabel('Sampling Instant'), ylabel('u_2(k)')
title('Model Validation: Manipulated Inputs');


%% MISO-II: Predictions on Validation data for ARMAX and OE Model 
% Infinite horizon predictions on validation data
n_st_ARMAX = length(phy2_ARMAX);
n_st_OE = length(phy2_OE);
xk_ARMAX = zeros(n_st_ARMAX, N_samples_val);
xk_OE = zeros(n_st_OE, N_samples_val);

% Predictions for validation for MISO-II ARMAX Model
yk2_hat_ARMAX(1) = C_mat2_ARMAX * xk_ARMAX(:,1);
for k = 1:N_samples_val-1
   xk_ARMAX(:,k+1) = phy2_ARMAX  * xk_ARMAX(:,k) + gama2_ARMAX * duk_val(:,k);
   yk2_hat_ARMAX(k+1) = C_mat2_ARMAX * xk_ARMAX(:,k+1);
end
yk2_hat_inf_ARMAX = yk2_hat_ARMAX;

% MISO-II: One Step Prediction on validation data ARMAX Model
xk_ARMAX = zeros(n_st_ARMAX, N_samples_val);
yk2_hat_ARMAX(1) = C_mat2_ARMAX * xk_ARMAX(:,1);
for k = 1: N_samples_val-1
    err_2(k) = dyk_val(2,k) - C_mat2_ARMAX * xk_ARMAX(:,k);
    xk_ARMAX(:,k+1) = phy2_ARMAX * xk_ARMAX(:,k) + gama2_ARMAX * duk_val(:,k) + K_mat2_ARMAX * err_2(k);
    yk2_hat_ARMAX(k+1) = C_mat2_ARMAX * xk_ARMAX(:,k+1) + err_2(k);
end
yk2_hat_one_ARMAX = yk2_hat_ARMAX;

% Predictions for validation for MISO-II OE Model
xk_OE = zeros(n_st_OE, N_samples_val);
yk2_hat_OE(1) = C_mat2_OE * xk_OE(:,1);
for k = 1: N_samples_val-1
    xk_OE(:,k+1) = phy2_OE * xk_OE(:,k) + gama2_OE * duk_val(:,k);
    yk2_hat_OE(k+1) = C_mat2_OE * xk_OE(:,k+1);
end

%% Graphical visualization and comaparision
figure, subplot(2,2,[1 2]), plot(dyk_val(2,:),'xr'),grid
hold on, plot(yk2_hat_inf_ARMAX,'-k'), grid
plot(yk2_hat_one_ARMAX,'-b'), grid
plot(yk2_hat_OE, 'g'), 
legend('Measure Output','Infinite Horizon Prediction','One Step Prediction','OE Prediction');
xlabel('Sampling Instant'), ylabel('Measured and Predicted Outputs')
title('Model Validation: MISO-II Predictions Comparisions')

subplot(223), stairs(duk_id(1,:)')
xlabel('Sampling Instant'),ylabel('u_1(k)')
title('Model Validation: Manipulated Inputs');
subplot(224), stairs(duk_id(2,:)')
xlabel('Sampling Instant'), ylabel('u_2(k)')
title('Model Validation: Manipulated Inputs');

%% ARMAX Model: Infinite horizon predictions on validation data for MIMO
n_st_ARMAX = length(phy_ARMAX);
n_st_OE = length(phy_OE);
xk_ARMAX = zeros(n_st_ARMAX, N_samples_val);
xk_OE = zeros(n_st_OE, N_samples_val);

% Predicttions for validation for MIMO ARMAX Model
yk_hat_ARMAX(:,1) = C_mat_ARMAX * xk_ARMAX(:,1);
for k = 1:N_samples_val-1
   xk_ARMAX(:,k+1) = phy_ARMAX  * xk_ARMAX(:,k) + gama_ARMAX * duk_val(:,k);
   yk_hat_ARMAX(:,k+1) = C_mat_ARMAX * xk_ARMAX(:,k+1);
end
yk_hat_inf_ARMAX = yk_hat_ARMAX;

% One Step Prediction on validation data for MIMO Model - ARMAX
xk_ARMAX = zeros(n_st_ARMAX, N_samples_val);
yk_hat_ARMAX(:,1) = C_mat_ARMAX * xk_ARMAX(:,1);
for k = 1: N_samples_val-1
    err(:,k) = dyk_val(:,k) - C_mat_ARMAX * xk_ARMAX(:,k);
    xk_ARMAX(:,k+1) = phy_ARMAX * xk_ARMAX(:,k) + gama_ARMAX * duk_val(:,k) + K_mat_ARMAX * err(:,k);
    yk_hat_ARMAX(:,k+1) = C_mat_ARMAX * xk_ARMAX(:,k+1) + err(:,k);
end
yk_hat_one_ARMAX = yk_hat_ARMAX;

% Predictions for validation for MIMO OE Model
xk_OE = zeros(n_st_OE, N_samples_val);
yk_hat_OE(:,1) = C_mat_OE * xk_OE(:,1);
for k = 1: N_samples_val-1
    xk_OE(:,k+1) = phy_OE * xk_OE(:,k) + gama_OE * duk_val(:,k);
    yk_hat_OE(:,k+1) = C_mat_OE * xk_OE(:,k+1);
end

%% Graphical Represenations for MIMO Model on Validation data
figure, subplot(211), plot(dyk_val(1,:),'xr'), grid, hold on
plot(yk_hat_inf_ARMAX(1,:), '-k'), grid
plot(yk_hat_one_ARMAX(1,:), '-b'), grid
plot(yk_hat_OE(1,:), '-g'), grid, hold off
legend('Measure Output','Infinite Horizon Prediction','One Step Prediction','OE Prediction');
xlabel('Sampling Instant'), ylabel('Measured and Predicted Outputs')
title('Model Validation: MIMO Predictions Comparisions of measurement  Y1')

subplot(212), plot(dyk_val(2,:),'xr'), grid, hold on
plot(yk_hat_inf_ARMAX(2,:), '-k'), grid
plot(yk_hat_one_ARMAX(2,:), '-b'), grid
plot(yk_hat_OE(2,:), '-g'), grid, hold off
legend('Measure Output','Infinite Horizon Prediction','One Step Prediction','OE Prediction');
xlabel('Sampling Instant'), ylabel('Measured and Predicted Outputs')
title('Model Validation: MIMO Predictions Comparisions of measurement Y2')

%% View step responses/ Nyquist plot using LTIVIEW
%ltiview(dss1_ARMAX, dss1_OE)
%ltiview(dss2_ARMAX, dss2_OE)
ltiview(dss_ARMAX, dss_OE)
