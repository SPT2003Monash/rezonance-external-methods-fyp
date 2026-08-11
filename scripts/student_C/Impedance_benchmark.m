clc; close all;clear all

%Obtaining the state space model. 
main_dir = 'C:\Users\erict\rezonance-external-methods-fyp\scripts\student_C';
cd(main_dir)
run("GFMI_Eric_Model.mlx")

% Run VSG_model.mlx to obtain all the individual state space matricies for
% individual componeents, then convert those matricies to transfer functions
main_dir='C:\Users\erict\rezonance-external-methods-fyp\scripts\student_C'; % set this to your parent directory-whereever you have this file
MatPower_path='C:\matpower8.1'; %Directry where you have Matpower installed

cd(main_dir)
addpath(genpath(MatPower_path));

GFMI_1 = Unified_GFMI;

load("Parameters_E1.mat")

%% PowerFlow steady state operating conditions. 
mpc = SMIB_PowerFlow(rg, lg);  
PowerFlow_results = runpf(mpc);

% Extract bus 1 angle (in degrees) and voltage magnitude (p.u). Angle
% relative to Grid
Angle = PowerFlow_results.bus(1, 9);         % Load angle in degrees
Voltage_mag = PowerFlow_results.bus(1, 8);    % Voltage magnitude in p.u.

% Calculate Line-to-Neutral Voltage Phasors in the Grid Reference Frame (Vg = 0 deg)
Angle_rad = deg2rad(Angle); % Positive leading angle
V_LN_RMS = (V_LL / sqrt(3));

Vc_phasor = Voltage_mag * V_LN_RMS * exp(1j * Angle_rad); % Vc leads Vg
Vg_phasor = V_LN_RMS * exp(1j * 0);                       % Grid aligned at 0 deg

% Calculate Complex Branch Currents in Grid Frame (RMS)
Zt = (rf2 + rg) + 1j * 2 * pi * 50 * (lf2 + lg);
Z2 = rf2 + 1j * 2 * pi * 50 * lf2;

I2_phasor = (Vc_phasor - Vg_phasor) / Zt; % Grid-side current phasor (RMS)
Vpcc_phasor = Vc_phasor - I2_phasor * Z2;

% Angle extracted from PSCAD. Found 1.75 to be the output, subtract pi/2
% due to cosine and sine swap? 
delta0  = 0; 
% Transform Phasors to Inverter d-q Frame (Aligned with Vc)
% To align frame with Vc, multiply by exp(-1j * Angle_rad) and sqrt(2) for peak
Vc_dq0 = Vc_phasor * exp(-1j * delta0) * sqrt(2); 
Vc_d0  = real(Vc_dq0); % Should equal Voltage_mag * V_LN_RMS * sqrt(2)
Vc_q0  = imag(Vc_dq0); % Should equal exactly 0.0

% Finding Initial operating conditions. 
I2_dq0 = I2_phasor * exp(-1j * delta0) * sqrt(2);
i2_d0  = real(I2_dq0); 
i2_q0  = imag(I2_dq0); 

% Store operating points

Vpcc_D0 = real(Vpcc_phasor * exp(-1j*0) * sqrt(2));
Vpcc_Q0 = imag(Vpcc_phasor * exp(-1j * 0) * sqrt(2));

% Validation for power calculations. 
Vpcc_d0 = real(Vpcc_phasor * exp(-1j*delta0) * sqrt(2));
Vpcc_q0 = imag(Vpcc_phasor * exp(-1j * delta0) * sqrt(2));

% Calculate Missing Steady-State Phasors (I1 and Vcf)
omega = 2 * pi * 50; % Nominal angular frequency

% Capacitor current (Ic) using Vc_phasor and filter capacitance (cf)
Ic_phasor = Vc_phasor * (1j * omega * cf);

% Inverter-side current (I1) using KCL: I1 = I2 + Ic
I1_phasor = I2_phasor + Ic_phasor;

% Inverter terminal voltage (Vcf) using KVL: Vcf = Vc + I1 * Z1
Z1 = rf1 + 1j * omega * lf1;
Vcf_phasor = Vc_phasor + I1_phasor * Z1;

% Transform Missing Phasors to Inverter d-q Frame
I1_dq0 = I1_phasor * exp(-1j * delta0) * sqrt(2);
i1_d0  = real(I1_dq0);
i1_q0  = imag(I1_dq0);

% Transform Vcf into inverter frame. 
Vcf_dq0 = Vcf_phasor * exp(-1j * delta0) * sqrt(2);
Vcf_d0  = real(Vcf_dq0);
Vcf_q0  = imag(Vcf_dq0);
%%
%operating point validation. Cross check these with PSCAD steady state to
realpowerInverter = 3/2*(Vc_d0*i2_d0 + Vc_q0*i2_q0)
reactiveInverter = 3/2*(Vc_q0*i2_d0-Vc_d0*i2_q0)

realpowerPCC = 3/2*(Vpcc_d0*i2_d0 + Vpcc_q0*i2_q0)
reactivePCC = 3/2*(Vpcc_q0*i2_d0-Vpcc_d0*i2_q0)
%% Closed Loop Substitution
% 1. List ALL symbolic targets for the FULL system
sym_targets = {'Vc_d0', 'Vc_q0', 'beta_i', 'beta_v', 'cf', 'cs', 'delta0', ...
               'i2_d0', 'i2_q0', 'k_ii', 'k_pi', 'k_iv', 'k_pv', ...
               'lf1', 'lf2', 'lg', 'rd', 'rf1', 'rf2', 'rg', 'rs', 'w', 'wn','td', ...
               'mp', 'wc_p', 'V_set', 'nq', 'wc', 'R_vir', 'X_vir', ... 
               'Vpcc_D0', 'Vpcc_Q0'}; 

% 2. List the exact workspace variables to substitute 
num_values  = {Vc_d0, Vc_q0, beta_i, beta_v, cf, cs, delta0, ...
               i2_d0, i2_q0, k_ii, k_pi, k_iv, k_pv, ...
               lf1, lf2, lg, rd, rf1, rf2, rg, rs, w, wn, td, ...
               mp, wc_p, V_set, nq, wc_q, R_vir, X_vir, ...
               Vpcc_D0, Vpcc_Q0}; 

% 3. Substitute values into the FULL Unified System Matrices

A_num = double(subs(A_sys0, sym_targets, num_values));
B_num = double(subs(B_sys0, sym_targets, num_values));
C_num = double(subs(C_sys0, sym_targets, num_values));
D_num = double(subs(D_sys0, sym_targets, num_values));

%% --- Create State-Space Model ---
% Build the full system
sys_full = ss(A_num, B_num, C_num, D_num);

%from the impedance state space model. 
% Inputs is U_sys = [p_ref; q_ref; Vpcc_D; Vpcc_Q];
% Outputs is Ystac = [i2_D; i2_Q; Vc_d; Vc_q; p_m; q_m];

% Admittance is 1/Z or I/V. Where I is the output and V is the input. 
%To get the matrices, just change the indices. 

% Extract Admittance Transfer Functions
% Input 3 = Vpcc_D, Input 4 = Vpcc_Q
% Output 1 = i2_D, Output 2 = i2_Q 
% Admittance Y = I_out / V_in
Y_dd = sys_full(1, 3); % Output 1 (i2_D) wrt Input 3 (Vpcc_D)
Y_qq = sys_full(2, 4); % Output 2 (i2_Q) wrt Input 4 (Vpcc_Q)
Y_dq = sys_full(1, 4); % Output 1 (i2_D) wrt Input 4 (Vpcc_Q)
Y_qd = sys_full(2, 3); % Output 2 (i2_Q) wrt Input 3 (Vpcc_D)
%% Reading CSV data
csvData = readtable('AVM_Full.csv', 'NumHeaderLines', 1);
% Extract frequencies and convert CSV magnitudes to dB
csv_f   = csvData.Var1;
csv_dd_db = 20*log10(abs(csvData.Var2));
csv_dq_db = 20*log10(abs(csvData.Var3));
csv_qd_db = 20*log10(abs(csvData.Var4));
csv_qq_db = 20*log10(abs(csvData.Var5));
%% --- Plotting the Admittance Matrix (Magnitudes Only) ---
frequencies_Hz = logspace(-3, 6, 100000); % Sweep from 1 Hz to 10000 Hz
frequencies_rad = 2 * pi * frequencies_Hz;

% 1. Get Bode magnitude data. (no phase) 
[mag_dd, ~] = bode(Y_dd, frequencies_rad);
[mag_dq, ~] = bode(Y_dq, frequencies_rad);
[mag_qd, ~] = bode(Y_qd, frequencies_rad);
[mag_qq, ~] = bode(Y_qq, frequencies_rad);

% 2. Convert magnitude to dB
mag_dd_db = 20*log10(squeeze(mag_dd*Zb)); 
mag_dq_db = 20*log10(squeeze(mag_dq*Zb));
mag_qd_db = 20*log10(squeeze(mag_qd*Zb));
mag_qq_db = 20*log10(squeeze(mag_qq*Zb));

% 3. Create the 2x2 Figure
figure('Color', 'w', 'Position', [100, 100, 1000, 700]); 
sgtitle('Analytical GFMI Admittance', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 16);


%% --- TOP LEFT: Y_dd Magnitude ---
subplot(2,2,1);
semilogx(frequencies_Hz, mag_dd_db, 'b', 'LineWidth', 1.5);
hold on
semilogx(csv_f, csv_dd_db); % CSV points
hold off
grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('d-d Axis Admittance (Y_{dd})');
xlim([1, 10000]);
ylim([-60, 60])

%% --- TOP RIGHT: Y_dq Magnitude ---
subplot(2,2,2);
semilogx(frequencies_Hz, mag_dq_db, 'g', 'LineWidth', 1.5);
hold on
semilogx(csv_f, csv_dq_db); % CSV points
hold off
grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('d-q Cross-Coupling Admittance (Y_{dq})');
xlim([1, 10000]);
ylim([-60, 60])
%% --- BOTTOM LEFT: Y_qd Magnitude ---
subplot(2,2,3);
semilogx(frequencies_Hz, mag_qd_db, 'm', 'LineWidth', 1.5);
hold on
semilogx(csv_f, csv_qd_db); % CSV points
hold off
grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('q-d Cross-Coupling Admittance (Y_{qd})');
xlim([1, 10000]);
ylim([-60, 60])
%% --- BOTTOM RIGHT: Y_qq Magnitude ---
subplot(2,2,4);
semilogx(frequencies_Hz, mag_qq_db, 'r', 'LineWidth', 1.5);
hold on
semilogx(csv_f, csv_qq_db); % CSV points
hold off
grid on; 
ylabel('Magnitude (dB)'); 
xlabel('Frequency (Hz)');
title('q-q Axis Admittance (Y_{qq})');
xlim([1, 10000]);
ylim([-60, 60])

%% Individual validation. LCL. filter. 
% 
% Vi_d_op = Vc_d0 + rf1*i1_d0 - lf1*i1_q0*w;
% Vi_q_op = Vc_q0 + rf1*i1_q0 + lf1*i1_d0*w;
% % 1. Define your operating point variables
% % (Replace these arbitrary numbers with your actual calculated operating points)
% op_vars = {'Vi_d', 'Vi_q', 'i1_d', 'i1_q', 'Vcf_d', 'Vcf_q', 'Vc_d', 'Vc_q', 'i2_d', 'i2_q', 'Vpcc_d', 'Vpcc_q', 'w_dev', 'lf1', 'rf1', 'rd', 'cf', 'lf2', 'rf2', 'w'};
% op_vals = [Vi_d_op,Vi_q_op,i1_d0,i1_q0,Vcf_d0,Vcf_q0,Vc_d0,Vc_q0,i2_d0,i2_q0,Vpcc_d0,Vpcc_q0,w_dev,lf1,rf1,rd,cf,lf2,rf2,w ]; % Your numeric values corresponding to the list above
% 
% % syms Vi_d Vi_q i1_d i1_q Vcf_d Vcf_q Vc_d Vc_q i2_d i2_q ig_d ig_q Vsypcc_d Vpcclf1_q w_dev %Variables
% % syms lf1 rf1 rd cf lf2 rf2 w %Parameters 
% % %Vectors
% % x_LCL = [i1_d; i1_q; Vcf_d; Vcf_q ; i2_d; i2_q]; 
% % e_LCL= [Vc_d; Vc_q]; 
% % u_LCL = [Vi_d; Vi_q; Vpcc_d; Vpcc_q;w_dev]; 
% % y_LCL = [i1_d; i1_q; Vc_d; Vc_q; i2_d; i2_q];
% 
% % 2. Substitute operating points to create numeric A, B, C, D matrices
% A_num = double(subs(A_LCL, op_vars, op_vals));
% B_num = double(subs(B_LCL, op_vars, op_vals));
% C_num = double(subs(C_LCL, op_vars, op_vals));
% D_num = double(subs(D_LCL, op_vars, op_vals));
% 
% % 3. Create the full Linear Time-Invariant (LTI) State-Space model
% sys_LCL = ss(A_num, B_num, C_num, D_num);
% 
% % 4. Extract the Admittance Matrix Y(s)
% % sys(outputs, inputs) -> sys([5, 6], [3, 4])
% Y_sys = -sys_LCL(5:6, 3:4); 
% 
% % 5. Configure Bode Plot Options (ensuring log scale and Hz if preferred)
% opts = bodeoptions('cstprefs');
% opts.FreqUnits = 'Hz'; 
% opts.MagUnits = 'dB';
% opts.PhaseUnits = 'deg';
% 
% % 6. Plot
% figure;
% bode(Y_sys, opts);
% grid on;
% title('LCL Filter Admittance Y(s)');
% 

