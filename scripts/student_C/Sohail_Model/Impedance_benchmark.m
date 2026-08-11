clc; close all; clear all

%Set directory and run the state space model file. 
main_dir = "C:\Users\erict\rezonance-external-methods-fyp\scripts\student_C\Sohail_Model";
cd(main_dir)
run("OL_Z_ConventionalPI.mlx")

%Create the state space model
GFMI = Unified_GFMI;

load("Parameters_Sohail.mat")

%% Find steady state conditions. 
%Analytical Steady-State Operating Point (P = 0, Q = 0)

% Grid voltage magnitude assumption
Vg_pu = 1.0; 

% 1. Filter output current (I2)
i2_d0 = 0;
i2_q0 = 0;

% 2. PCC Voltages (Aligned with Grid)
Vpcc_d0 = Vg_pu;
Vpcc_q0 = 0;

% 3. Capacitor Voltages (accounting for damping resistor rd)
w_Cf_Rd = cf_pu * rd_pu;

vc_d0 = Vpcc_d0 / (1 + w_Cf_Rd^2);
vc_q0 = (-Vpcc_d0 * w_Cf_Rd) / (1 + w_Cf_Rd^2);

% 4. Misalignment angle
Delta0 = 0; 

% Duplicate variables if your symbolic model specifically calls for capital letters
i2_D0 = i2_d0;
i2_Q0 = i2_q0;

V_nom =1;
%% Sub in numerical values. 
% 1. List ALL symbolic targets for the FULL system
sym_targets = {'lf1_pu', 'rf1_pu', 'cf_pu', 'rd_pu', 'td', ...
               'k_pi_pu', 'k_ii_pu', 'k_pv_pu', 'k_iv_pu', 'beta_v', 'beta_i', ...
               'KpQ_pu', 'KiQ_pu', 'mp', 'H', 'V_nom', ...
               'w_base', 'w', ...
               'Vc_d0', 'Vc_q0', 'i2_D0', 'i2_Q0', 'i2_d0', 'i2_q0', ...
               'rv_pu', 'lv_pu'};
% 2. List the exact workspace variables to substitute 
num_values  = {lf1_pu, rf1_pu, cf_pu,rd_pu,td,...
               k_pi_pu,k_ii_pu,k_pv_pu,k_iv_pu,beta_v,beta_i,...
               KpQ_pu,KiQ_pu,mp,H,V_nom,...
               w_base,w,...
               vc_d0,vc_q0,i2_D0,i2_Q0,i2_d0,i2_q0,...
               rv_pu,lv_pu}; 

% 3. Substitute values into the FULL Unified System Matrices

A_num = double(subs(A_sys0, sym_targets, num_values));
B_num = double(subs(B_sys0, sym_targets, num_values));
C_num = double(subs(C_sys0, sym_targets, num_values));
D_num = double(subs(D_sys0, sym_targets, num_values));

%% Create admittance matrix. 
Y_sys = ss(A_num, B_num, C_num, D_num);

Y_dd = Y_sys(1, 1);
Y_qq = Y_sys(2, 2);
Y_dq = Y_sys(1, 2);
Y_qd = Y_sys(2, 1);


%% --- Plotting the Admittance Matrix (Magnitudes Only) ---
frequencies_Hz = logspace(-3, 6, 100000); % Sweep from 1 Hz to 10000 Hz
frequencies_rad = 2 * pi * frequencies_Hz;

% 1. Get Bode magnitude data. (no phase) 
[mag_dd, ~] = bode(Y_dd, frequencies_rad);
[mag_dq, ~] = bode(Y_dq, frequencies_rad);
[mag_qd, ~] = bode(Y_qd, frequencies_rad);
[mag_qq, ~] = bode(Y_qq, frequencies_rad);

% 2. Convert magnitude to dB
mag_dd_db = 20*log10(squeeze(mag_dd)); 
mag_dq_db = 20*log10(squeeze(mag_dq));
mag_qd_db = 20*log10(squeeze(mag_qd));
mag_qq_db = 20*log10(squeeze(mag_qq));

% 3. Create the 2x2 Figure
figure('Color', 'w', 'Position', [100, 100, 1000, 700]); 
sgtitle('Analytical GFMI Admittance', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 16);


%% --- TOP LEFT: Y_dd Magnitude ---
subplot(2,2,1);
semilogx(frequencies_Hz, mag_dd_db, 'b', 'LineWidth', 1.5);

grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('d-d Axis Admittance (Y_{dd})');
xlim([1, 10000]);
ylim([-60, 60])

%% --- TOP RIGHT: Y_dq Magnitude ---
subplot(2,2,2);
semilogx(frequencies_Hz, mag_dq_db, 'g', 'LineWidth', 1.5);

grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('d-q Cross-Coupling Admittance (Y_{dq})');
xlim([1, 10000]);
ylim([-60, 60])
%% --- BOTTOM LEFT: Y_qd Magnitude ---
subplot(2,2,3);
semilogx(frequencies_Hz, mag_qd_db, 'm', 'LineWidth', 1.5);

grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('q-d Cross-Coupling Admittance (Y_{qd})');
xlim([1, 10000]);
ylim([-60, 60])
%% --- BOTTOM RIGHT: Y_qq Magnitude ---
subplot(2,2,4);
semilogx(frequencies_Hz, mag_qq_db, 'r', 'LineWidth', 1.5);

grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('q-q Axis Admittance (Y_{qq})');
xlim([1, 10000]);
ylim([-60, 60])