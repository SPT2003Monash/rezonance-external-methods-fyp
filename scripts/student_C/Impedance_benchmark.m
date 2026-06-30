clc; close all;clear all

main_dir = 'C:\Users\erict\rezonance-external-methods-fyp\scripts\student_C';
cd(main_dir)
run("GFMI_Eric_Model.mlx")

%Run VSG_model.mlx to obtain all the individual state space matricies for
%individual componeents, then convert those matricies to transfer functions



main_dir='C:\Users\erict\rezonance-external-methods-fyp\scripts\student_C'; % set this to your parent directory-whereever you have this file
MatPower_path='C:\matpower8.1'; %Directry where you have Matpower installed

cd(main_dir)
addpath(genpath(MatPower_path));

GFMI_1 = Unified_GFMI;

load("Parameters.mat")





%%
%run powerflow analysis for Vdq0 and Idq0 and angle0
mpc = SMIB_PowerFlow(rg,lg);  

% Run the power flow for this operating point
PowerFlow_results = runpf(mpc);

% Extract bus 1 angle (in degrees) and voltage magnitude (p.u)
Angle = PowerFlow_results.bus(1,9);
Voltage_mag = PowerFlow_results.bus(1,8);

% Calculate space phasors
Angle_rad = deg2rad(Angle);
Vc_phasor = Voltage_mag * V_LL/sqrt(3) * exp(1j * Angle_rad); % line-to-neutral voltage
Vg_phasor = V_LL/sqrt(3);

% Compute impedance values (assuming these parameters are constant)
Zt = (rf2 + rg + 1j * 2*pi*50 * (lf2 + lg));
Z2 = rf2 + 1j * 2*pi*50 * lf2;

% Calculate current and PCC voltage phasors
I2_phasor = (Vc_phasor - Vg_phasor) / Zt;
Vpcc_phasor = (Vc_phasor - I2_phasor * Z2) * sqrt(2);

% Store operating point values
delta0 = Angle_rad;    
Vpcc_D0 = real(Vpcc_phasor);
Vpcc_Q0 = imag(Vpcc_phasor);

% Transform Vc phasor into d–q components
Vc_dq0 = Vc_phasor * exp(-1j * Angle_rad) * sqrt(2);
Vc_d0 = real(Vc_dq0);
Vc_q0 = imag(Vc_dq0);

% Transform I2 phasor into d–q components
I2_dq0 = I2_phasor * exp(-1j * Angle_rad) * sqrt(2);
i2_d0 = real(I2_dq0);
i2_q0 = imag(I2_dq0);

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
% Inputs is U_sys = [p_ref; q_ref; Vg_D; Vg_Q];
% Outputs is Ystac = [i1_d; i1_q; Vc_d; Vc_q; p_m; q_m];

% Admittance is 1/Z or I/V. Where I is the output and V is the input. 
%To get the matrices, just change the indices. 


% Extract Admittance Transfer Functions
% Input 3 = Vg_D, Input 4 = Vg_Q
% Output 1 = i1_d, Output 2 = i1_q 
Y_dd = sys_full(1, 3); % D-axis Admittance (Delta i1_d / Delta Vg_D)
Y_qq = sys_full(2, 4); % Q-axis Admittance (Delta i1_q / Delta Vg_Q)
Y_dq = sys_full(1, 4); % Cross-coupling (Delta i1_d / Delta Vg_Q)
Y_qd = sys_full(2, 3); % Cross-coupling (Delta i1_q / Delta Vg_D)

%% --- Plotting the Admittance Matrix (Magnitudes Only) ---
frequencies_Hz = logspace(0, 4, 1000); % Sweep from 1 Hz to 10000 Hz
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










%Extended analysis for later parts of the project. 
% %% --- Eigenvalue Stability Analysis (Pole Map) ---
% % 1. Calculate the Eigenvalues of the closed-loop system
% % The A_num matrix perfectly dictates absolute stability
% sys_eig = eig(A_num);
% 
% % 2. Create the Pole Map Figure
% figure('Color', 'w', 'Position', [200, 200, 800, 600]);
% scatter(real(sys_eig), imag(sys_eig)/(2*pi), 120, 'bx', 'LineWidth', 2);
% hold on; grid on;
% 
% % 3. Draw Stability Boundaries
% % The Y-axis is the boundary between life and death (Stable vs Unstable)
% xline(0, 'r--', 'LineWidth', 2, 'DisplayName', 'Instability Boundary'); 
% yline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
% 
% % 4. Formatting the Plot
% xlabel('Real Part / Damping (\sigma)');
% ylabel('Imaginary Part / Oscillation Frequency (Hz)');
% title('Closed-Loop System Eigenvalues (Stability Pole Map)');
% xlim([-2000, 200]); % Focused zoom on the dominant low-frequency poles
% legend('System Poles', 'Location', 'northeast');
% 
% % 5. Mathematical Stability Report (Prints to Command Window)
% fprintf('\n========================================\n');
% fprintf('       SYSTEM STABILITY REPORT\n');
% fprintf('========================================\n');
% 
% % Ignore poles exactly at 0 (integrators/reference frame alignment)
% non_zero_eigs = sys_eig(abs(real(sys_eig)) > 1e-5);
% 
% % Find the pole closest to the Right-Half Plane (The Dominant Pole)
% [~, min_idx] = max(real(non_zero_eigs)); 
% dominant_pole = non_zero_eigs(min_idx);
% 
% % Check for Right-Half Plane (RHP) poles
% if real(dominant_pole) > 0
%     fprintf('STATUS: UNSTABLE (Poles detected in Right-Half Plane!)\n');
% else
%     fprintf('STATUS: STABLE (All dynamic poles are in the Left-Half Plane)\n');
% end
% 
% % Calculate specific metrics for the dominant mode
% oscillation_freq_Hz = abs(imag(dominant_pole)) / (2*pi);
% damping_ratio = -real(dominant_pole) / abs(dominant_pole);
% 
% fprintf('Dominant Pole Location : %.2f +/- %.2fj Hz\n', real(dominant_pole), oscillation_freq_Hz);
% fprintf('Damping Ratio (\zeta)    : %.4f\n', damping_ratio);
% fprintf('========================================\n\n');
% 
% 
% %% --- Open-Loop Determinant Analysis (det(Y)) ---
% % We use freqresp to get the complex numerical values at each frequency 
% % to avoid massive state-space inflation.
% Y_dd_resp = squeeze(freqresp(Y_dd, frequencies_rad));
% Y_qq_resp = squeeze(freqresp(Y_qq, frequencies_rad));
% Y_dq_resp = squeeze(freqresp(Y_dq, frequencies_rad));
% Y_qd_resp = squeeze(freqresp(Y_qd, frequencies_rad));
% 
% % Initialize the determinant array
% det_Y = zeros(size(frequencies_rad));
% 
% % Calculate the determinant at each frequency step
% % det(Y) = (Y_dd * Y_qq) - (Y_dq * Y_qd)
% for k = 1:length(frequencies_rad)
%     det_Y(k) = Y_dd_resp(k) * Y_qq_resp(k) - Y_dq_resp(k) * Y_qd_resp(k);
% end
% 
% % Convert the determinant to Magnitude (dB) and Phase (degrees)
% mag_det_db = 20*log10(abs(det_Y));
% % unwrap() prevents the phase from jumping abruptly from +180 to -180
% phase_det_deg = unwrap(angle(det_Y)) * (180/pi); 
% 
% %% --- Plot 1: Bode Plot of the Determinant ---
% figure('Color', 'w', 'Position', [250, 150, 800, 600]);
% 
% % Magnitude Plot
% subplot(2,1,1);
% semilogx(frequencies_Hz, mag_det_db, 'k', 'LineWidth', 1.5);
% grid on;
% ylabel('Magnitude (dB)');
% title('Determinant of Inverter Admittance Matrix: \det(Y)');
% xlim([1, 10000]);
% 
% % Phase Plot
% subplot(2,1,2);
% semilogx(frequencies_Hz, phase_det_deg, 'k', 'LineWidth', 1.5);
% grid on;
% ylabel('Phase (deg)');
% xlabel('Frequency (Hz)');
% xlim([1, 10000]);
% 
% %% --- Plot 2: Nyquist Plot of the Determinant ---
% % Often used for the Generalized Nyquist Criterion (GNC) check
% figure('Color', 'w', 'Position', [300, 200, 600, 600]);
% plot(real(det_Y), imag(det_Y), 'b', 'LineWidth', 1.5); hold on;
% % Plot the negative frequency conjugate mirroring for a complete Nyquist contour
% plot(real(det_Y), -imag(det_Y), 'b--', 'LineWidth', 1.2); 
% plot(0, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Origin reference point
% grid on;
% 
% % Formatting
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('Nyquist Plot of Admittance Determinant');
% xline(0, 'k-', 'HandleVisibility', 'off'); 
% yline(0, 'k-', 'HandleVisibility', 'off');
% legend('Positive Frequencies', 'Negative Frequencies', 'Location', 'best');
% 
% 
% %% --- Time-Domain Simulation Setup ---
% % 1. Create the continuous-time state-space object
% sys_cl = ss(A_num, B_num, C_num, D_num);
% 
% % 2. Define the Time Vector
% t_end = 0.5;                  % Simulate for 0.5 seconds
% dt = 1/f_step;                % Match your parameter step time
% t = 0:dt:t_end;
% N = length(t);
% 
% % 3. Define the Input Matrix (U_sys = [p_ref; q_ref; Vg_D; Vg_Q])
% % Initialize all inputs to zero deviation (starts exactly at your P_inj = 0 operating point)
% u_sim = zeros(N, 4); 
% 
% % --- TEST CASE: 0 to 1.0 pu Step Increase in Active Power Reference ---
% % Since baseline P_inj = 0, the initial deviation is 0.
% % At t = 0.1 seconds, we apply a full 1.0 pu step perturbation.
% step_time = 0.1;
% step_idx = t >= step_time;
% u_sim(step_idx, 1) = 0.1 * S_base; % Step perturbation of 1.0 pu Watts
% 
% %% --- Run the Linear Simulation ---
% % lsim simulates the small-signal perturbations (Delta y)
% [delta_y, t_out, delta_x] = lsim(sys_cl, u_sim, t);
% 
% %% --- Reconstruct Total Physical Values (Large-Signal) ---
% % Reconstruct total physical waveforms by adding the baseline operating points
% i1_d_total  = i2_d0  + delta_y(:, 1); % Total D-axis Inverter Current
% i1_q_total  = i2_q0  + delta_y(:, 2); % Total Q-axis Inverter Current
% Vc_d_total  = Vc_d0  + delta_y(:, 3); % Total D-axis Cap Voltage
% Vc_q_total  = Vc_q0  + delta_y(:, 4); % Total Q-axis Cap Voltage
% P_meas_total = P_inj*S_base + delta_y(:, 5); % Total Measured Active Power (P_inj should be 0)
% Q_meas_total = Q_inj*S_base + delta_y(:, 6); % Total Measured Reactive Power
% 
% %% --- Plotting the Transient Response ---
% figure('Color', 'w', 'Position', [100, 100, 900, 700]);
% 
% % Plot Active Power Tracking
% subplot(3,1,1);
% plot(t, P_meas_total / S_base, 'b', 'LineWidth', 2); hold on;
% plot(t, (P_inj*S_base + u_sim(:,1)) / S_base, 'r--', 'LineWidth', 1.5);
% grid on; ylabel('Active Power (pu)');
% title('Time-Domain Step Response (0.9 to 1.0 pu P_{ref} Step Change)');
% legend('Measured P', 'Reference P', 'Location', 'southeast');
% ylim([-0.1, 1.2]); % Adjusted bounds to fit the full 0-1 range cleanly
% 
% % Plot Inverter Currents
% subplot(3,1,2);
% plot(t, i1_d_total, 'b', 'LineWidth', 1.5); hold on;
% plot(t, i1_q_total, 'r', 'LineWidth', 1.5);
% grid on; ylabel('Currents (A)');
% legend('i_{1d} (Active)', 'i_{1q} (Reactive)');
% 
% % Plot Capacitor Voltages
% subplot(3,1,3);
% plot(t, Vc_d_total, 'b', 'LineWidth', 1.5); hold on;
% plot(t, Vc_q_total, 'r', 'LineWidth', 1.5);
% grid on; ylabel('Voltages (V)'); xlabel('Time (seconds)');
% legend('V_{cd}', 'V_{cq}');