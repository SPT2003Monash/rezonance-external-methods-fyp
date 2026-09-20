% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Examines the effect of measurement-window length on PSCAD RL impedance and Richardson estimates.

clear; clc; close all;

%% Load PSCAD CSV data
% CSV files should have:
% Column 1 = time
% Column 2 = signal value

current_data = readmatrix('Ia.csv');
voltage_data = readmatrix('Ea.csv');

% Remove rows with NaN values, for example header rows
current_data = current_data(all(~isnan(current_data(:,1:2)),2),:);
voltage_data = voltage_data(all(~isnan(voltage_data(:,1:2)),2),:);

% Extract time and signal values
t_i = current_data(:,1);
i_raw = current_data(:,2);

t_v = voltage_data(:,1);
v_raw = voltage_data(:,2);

% Use current time vector as the main time vector
t = t_i;
i = i_raw;

% If voltage time vector is slightly different, interpolate voltage onto current time vector
if length(t_v) ~= length(t_i) || max(abs(t_v - t_i)) > 1e-9
    v = interp1(t_v, v_raw, t, 'linear', 'extrap');
else
    v = v_raw;
end
% Correct voltage polarity from PSCAD measurement direction
v = -v;
%% add artificial noise to PSCAD waveforms
rng(1);  % repeatable noise

noise_level = 0.001;  % 0.1% noise

i = i + noise_level * max(abs(i)) * randn(size(i));
v = v + noise_level * max(abs(v)) * randn(size(v));
%% PSCAD RL circuit parameters
% These must match the PSCAD circuit values

R = 1;          % Resistance in ohms
L = 0.1;        % Inductance in H
f_inj = 1;      % Injection frequency in Hz
w = 2*pi*f_inj;

% Analytical true impedance of RL circuit
Z_true = R + 1j*w*L;

%% Plot raw PSCAD waveforms
figure;
plot(t, i, 'LineWidth', 1.2); hold on;
plot(t, v, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Signal value');
legend('Current Ia (A)', 'Voltage Ea (V)', 'Location', 'best');
title('PSCAD RL Circuit Time-Domain Output');

%% scan-window settings
% Here, short/medium/long windows are created from the same PSCAD run.
% A shorter window uses fewer cycles, so it imitates a faster but lower-quality scan.
% A longer window uses more cycles, so it imitates a slower but higher-quality scan.

N_short  = 2;    % Short-window estimate uses last 2 cycles
N_medium = 5;    % Medium-window estimate uses last 5 cycles
N_long   = 8;    % Long-window estimate uses last 8 cycles

% h is treated as inverse of number of cycles.
% More cycles = smaller h = better scan estimate.
h_short  = 1/N_short;
h_medium = 1/N_medium;
h_long   = 1/N_long;

p = 2;           % Assumed convergence order

% Refinement ratios
r_SM = h_short / h_medium;   % short to medium
r_ML = h_medium / h_long;    % medium to long
r_SL = h_short / h_long;     % short to long

%% Extract impedance estimates from PSCAD data
% Each estimate comes from a different time-window length.

Z_short  = extract_impedance_from_window(t, v, i, f_inj, N_short);
Z_medium = extract_impedance_from_window(t, v, i, f_inj, N_medium);
Z_long   = extract_impedance_from_window(t, v, i, f_inj, N_long);

%% Apply Richardson refinement
% Richardson combines two estimates to try to cancel systematic error.

Z_R_SM = (r_SM^p * Z_medium - Z_short) / (r_SM^p - 1);
Z_R_ML = (r_ML^p * Z_long   - Z_medium) / (r_ML^p - 1);
Z_R_SL = (r_SL^p * Z_long   - Z_short) / (r_SL^p - 1);

%% Calculate relative errors
% Since this is a simple RL circuit, the analytical impedance is known.

err_short  = abs(Z_short  - Z_true) / abs(Z_true);
err_medium = abs(Z_medium - Z_true) / abs(Z_true);
err_long   = abs(Z_long   - Z_true) / abs(Z_true);

err_R_SM = abs(Z_R_SM - Z_true) / abs(Z_true);
err_R_ML = abs(Z_R_ML - Z_true) / abs(Z_true);
err_R_SL = abs(Z_R_SL - Z_true) / abs(Z_true);

%% Simple Richardson acceptance / rejection rule
% Accept Richardson only if it improves compared with the better input estimate.

accept_R_SM = err_R_SM < err_medium;
accept_R_ML = err_R_ML < err_long;
accept_R_SL = err_R_SL < err_long;

%% Plot relative error comparison
method_names = categorical({ ...
    'Short scan', ...
    'Medium scan', ...
    'Long scan', ...
    'Richardson S+M', ...
    'Richardson M+L', ...
    'Richardson S+L'});

method_names = reordercats(method_names, { ...
    'Short scan', ...
    'Medium scan', ...
    'Long scan', ...
    'Richardson S+M', ...
    'Richardson M+L', ...
    'Richardson S+L'});

errors = [err_short, err_medium, err_long, err_R_SM, err_R_ML, err_R_SL];

figure;
bar(method_names, errors);
grid on;
ylabel('Relative Error');
title('Week 10 PSCAD RL Test: Richardson Window Error Comparison');

%% Print summary table
fprintf('\n================ Week 10 PSCAD RL Richardson Window Test Summary ================\n');

fprintf('\nAnalytical true impedance:\n');
fprintf('Z_true = %.6f + j%.6f ohm\n', real(Z_true), imag(Z_true));
fprintf('|Z_true| = %.6f ohm\n', abs(Z_true));

fprintf('\nEstimated impedance values:\n');
fprintf('Short-window estimate:             %.6f + j%.6f ohm\n', real(Z_short),  imag(Z_short));
fprintf('Medium-window estimate:            %.6f + j%.6f ohm\n', real(Z_medium), imag(Z_medium));
fprintf('Long-window estimate:              %.6f + j%.6f ohm\n', real(Z_long),   imag(Z_long));
fprintf('Richardson short + medium:         %.6f + j%.6f ohm\n', real(Z_R_SM),   imag(Z_R_SM));
fprintf('Richardson medium + long:          %.6f + j%.6f ohm\n', real(Z_R_ML),   imag(Z_R_ML));
fprintf('Richardson short + long:           %.6f + j%.6f ohm\n', real(Z_R_SL),   imag(Z_R_SL));

fprintf('\nRelative Error:\n');
fprintf('Short-window scan:                 %.4e\n', err_short);
fprintf('Medium-window scan:                %.4e\n', err_medium);
fprintf('Long-window scan:                  %.4e\n', err_long);
fprintf('Richardson short + medium:         %.4e\n', err_R_SM);
fprintf('Richardson medium + long:          %.4e\n', err_R_ML);
fprintf('Richardson short + long:           %.4e\n', err_R_SL);

fprintf('\nRichardson Acceptance Check:\n');
fprintf('Short + medium accepted:           %d\n', accept_R_SM);
fprintf('Medium + long accepted:            %d\n', accept_R_ML);
fprintf('Short + long accepted:             %d\n', accept_R_SL);

fprintf('\nConclusion:\n');

if err_R_SM < err_medium
    fprintf('Richardson short + medium improved compared with the medium-window estimate.\n');
else
    fprintf('Richardson short + medium did NOT improve compared with the medium-window estimate.\n');
end

if err_R_ML < err_long
    fprintf('Richardson medium + long improved compared with the long-window estimate.\n');
else
    fprintf('Richardson medium + long did NOT improve compared with the long-window estimate.\n');
end

if err_R_SL < err_long
    fprintf('Richardson short + long improved compared with the long-window estimate.\n');
else
    fprintf('Richardson short + long did NOT improve compared with the long-window estimate.\n');
end

fprintf('This PSCAD RL test checks whether the synthetic Richardson workflow also works on real PSCAD time-domain data.\n');
fprintf('================================================================================\n');

%% Sign convention warning
% If the imaginary part is negative instead of positive, the current or voltage
% measurement direction may be reversed in PSCAD.

if imag(Z_long) < 0
    fprintf('\nWARNING: Imaginary part of Z_long is negative.\n');
    fprintf('For an RL circuit, expected impedance is R + jwL with positive imaginary part.\n');
    fprintf('This may mean the voltage polarity or current direction is reversed in PSCAD.\n');
end

%% Local function: impedance extraction from time window
function Z_est = extract_impedance_from_window(t, v, i, f_inj, N_cycles)

    % Angular frequency
    w = 2*pi*f_inj;

    % Period of injected sinusoid
    T = 1/f_inj;

    % Use the last N_cycles of the simulation
    t_end = max(t);
    t_start = t_end - N_cycles*T;

    idx = t >= t_start;

    t_win = t(idx);
    v_win = v(idx);
    i_win = i(idx);

    % Remove DC offsets
    v_win = v_win - mean(v_win);
    i_win = i_win - mean(i_win);

    % Reference sine and cosine waves
    cos_ref = cos(w*t_win);
    sin_ref = sin(w*t_win);

    % Extract current phasor
    I_cos = 2/length(t_win) * sum(i_win .* cos_ref);
    I_sin = 2/length(t_win) * sum(i_win .* sin_ref);

    % Extract voltage phasor
    V_cos = 2/length(t_win) * sum(v_win .* cos_ref);
    V_sin = 2/length(t_win) * sum(v_win .* sin_ref);

    % Complex phasors
    I_phasor = I_cos - 1j*I_sin;
    V_phasor = V_cos - 1j*V_sin;

    % Impedance estimate
    Z_est = V_phasor / I_phasor;
end