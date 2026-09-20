% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Validates pulse-spectrum acceptance criteria using an analytical RL circuit.

clear; clc; close all;

%% R-L benchmark parameters
R = 1;          % ohm
L = 0.1;        % H

%% Simulation settings
Tsim = 20;          % seconds
dt = 50e-6;         % sampling/logging step
fs = 1/dt;          % sampling frequency
t = (0:dt:Tsim-dt).';
N = length(t);

%% Pulse settings
pulse_start = 2;        % seconds
pulse_width = 0.2;      % seconds
pulse_amp = 0.01;       % A

%% Noise setting
noise_level = 0.001;    % 0.001 = 0.1% noise
rng(1);

%% Acceptance rule
threshold_ratio = 0.01;     % 1% of max input FFT magnitude
fmax_plot = 50;             % Hz

%% Frequencies for selected result table
% Low-frequency points and few extra screening points
selected_freqs = [0.2 0.5 1 2 5 10 20];

%% Pulse shapes to compare
pulse_shapes = ["rectangular", "raised_cosine", "gaussian"];

%% Frequency vector using signed frequency for analytical impedance
f_full = (0:N-1).' * fs/N;
f_signed = f_full;
f_signed(f_signed > fs/2) = f_signed(f_signed > fs/2) - fs;

Z_full = R + 1j*2*pi*f_signed*L;

%% Positive frequency vector
half_N = floor(N/2);
f_pos = f_full(1:half_N);
idx_plot = (f_pos > 0) & (f_pos <= fmax_plot);

%% Storage for pulse-shape summary
Shape = strings(length(pulse_shapes),1);
Accepted_Selected = zeros(length(pulse_shapes),1);
Rejected_Selected = zeros(length(pulse_shapes),1);
Accepted_0_to_50Hz = zeros(length(pulse_shapes),1);
Rejected_0_to_50Hz = zeros(length(pulse_shapes),1);
Mean_Complex_Error_Accepted_Selected = zeros(length(pulse_shapes),1);
Best_Selected_Frequency_Hz = zeros(length(pulse_shapes),1);
Best_Selected_Complex_Error_percent = zeros(length(pulse_shapes),1);

%% Main loop over pulse shapes
for s = 1:length(pulse_shapes)

    shape = pulse_shapes(s);
    Shape(s) = shape;

    fprintf('\n=============================================\n');
    fprintf('Processing pulse shape: %s\n', shape);

    %% Create pulse
    i_time = zeros(size(t));

    idx_pulse = (t >= pulse_start) & (t <= pulse_start + pulse_width);
    tau = (t(idx_pulse) - pulse_start) / pulse_width;

    switch shape

        case "rectangular"
            i_time(idx_pulse) = pulse_amp;

        case "raised_cosine"
            % Smooth pulse that starts and ends at zero
            i_time(idx_pulse) = pulse_amp * 0.5 .* (1 - cos(2*pi*tau));

        case "gaussian"
            % Gaussian-like finite pulse inside the same pulse window
            centre = pulse_start + pulse_width/2;
            sigma = pulse_width/6;
            i_time(idx_pulse) = pulse_amp * exp(-0.5*((t(idx_pulse)-centre)/sigma).^2);

        otherwise
            error('Unknown pulse shape.');
    end

    %% FFT of input current pulse
    I_fft = fft(i_time);

    %% Analytical voltage response in frequency domain
    % V(f) = Z(f) I(f)
    V_fft = Z_full .* I_fft;

    %% Convert voltage to time domain
    v_time = real(ifft(V_fft));

    %% Add measurement noise
    v_time_noisy = v_time + noise_level * max(abs(v_time)) * randn(size(v_time));
    i_time_noisy = i_time + noise_level * max(abs(i_time)) * randn(size(i_time));

    %% FFT of measured noisy signals
    V_meas_fft = fft(v_time_noisy);
    I_meas_fft = fft(i_time_noisy);

    V_pos = V_meas_fft(1:half_N);
    I_pos = I_meas_fft(1:half_N);

    %% True impedance and estimated impedance at positive frequencies
    Z_true_pos = R + 1j*2*pi*f_pos*L;
    Z_est_pos = V_pos ./ I_pos;

    %% Accept/reject based on input spectral energy
    I_mag = abs(I_pos);
    I_mag_no_dc = I_mag;
    I_mag_no_dc(1) = 0;

    threshold = threshold_ratio * max(I_mag_no_dc);

    accepted = I_mag > threshold;
    rejected = ~accepted;

    % Always reject DC
    accepted(1) = false;
    rejected(1) = true;

    %% Complex error at all positive bins
    complex_error_all = abs(Z_est_pos - Z_true_pos) ./ abs(Z_true_pos) * 100;

    %% Selected frequency table for this pulse shape
    selected_bin_freqs = zeros(length(selected_freqs),1);
    Ztrue_selected = zeros(length(selected_freqs),1);
    Zest_selected = zeros(length(selected_freqs),1);
    input_strength_selected = zeros(length(selected_freqs),1);
    mag_err_selected = zeros(length(selected_freqs),1);
    phase_err_selected = zeros(length(selected_freqs),1);
    complex_err_selected = zeros(length(selected_freqs),1);
    decision_selected = strings(length(selected_freqs),1);

    for k = 1:length(selected_freqs)

        f_target = selected_freqs(k);

        [~, idx_near] = min(abs(f_pos - f_target));

        selected_bin_freqs(k) = f_pos(idx_near);
        Ztrue_selected(k) = Z_true_pos(idx_near);
        Zest_selected(k) = Z_est_pos(idx_near);

        input_strength_selected(k) = I_mag(idx_near) / max(I_mag_no_dc) * 100;

        mag_err_selected(k) = abs(abs(Zest_selected(k)) - abs(Ztrue_selected(k))) ...
            / abs(Ztrue_selected(k)) * 100;

        phase_est = angle(Zest_selected(k))*180/pi;
        phase_true = angle(Ztrue_selected(k))*180/pi;
        phase_err_selected(k) = abs(wrapTo180_local(phase_est - phase_true));

        complex_err_selected(k) = abs(Zest_selected(k) - Ztrue_selected(k)) ...
            / abs(Ztrue_selected(k)) * 100;

        if accepted(idx_near)
            decision_selected(k) = "Accepted";
        else
            decision_selected(k) = "Rejected";
        end
    end

    PulseTable = table(selected_freqs.', selected_bin_freqs, ...
        real(Ztrue_selected), imag(Ztrue_selected), ...
        real(Zest_selected), imag(Zest_selected), ...
        input_strength_selected, ...
        mag_err_selected, phase_err_selected, complex_err_selected, ...
        decision_selected, ...
        'VariableNames', {'Target_Frequency_Hz', 'Nearest_FFT_Bin_Hz', ...
        'Ztrue_real', 'Ztrue_imag', ...
        'Zest_real', 'Zest_imag', ...
        'Input_FFT_strength_percent', ...
        'Magnitude_error_percent', 'Phase_error_deg', 'Complex_error_percent', ...
        'Decision'});

    fprintf('\nSelected frequency table for %s pulse:\n', shape);
    disp(PulseTable);

    %% Store summary values
    accepted_selected = decision_selected == "Accepted";
    rejected_selected = decision_selected == "Rejected";

    Accepted_Selected(s) = sum(accepted_selected);
    Rejected_Selected(s) = sum(rejected_selected);

    Accepted_0_to_50Hz(s) = sum(accepted(idx_plot));
    Rejected_0_to_50Hz(s) = sum(rejected(idx_plot));

    if any(accepted_selected)
        Mean_Complex_Error_Accepted_Selected(s) = mean(complex_err_selected(accepted_selected));
        [best_err, best_idx_local] = min(complex_err_selected(accepted_selected));
        accepted_freqs_temp = selected_freqs(accepted_selected);
        Best_Selected_Frequency_Hz(s) = accepted_freqs_temp(best_idx_local);
        Best_Selected_Complex_Error_percent(s) = best_err;
    else
        Mean_Complex_Error_Accepted_Selected(s) = NaN;
        Best_Selected_Frequency_Hz(s) = NaN;
        Best_Selected_Complex_Error_percent(s) = NaN;
    end

    %% Clean shape name for file names
    shape_file = char(shape);

    %% Plot 1: pulse in time domain
    figure;
    plot(t, i_time, 'LineWidth', 1.2);
    grid on;
    xlabel('Time (s)');
    ylabel('Current pulse i(t) (A)');
    title(sprintf('%s current pulse', strrep(shape_file, '_', '-')));
    xlim([pulse_start-0.5, pulse_start+pulse_width+0.5]);
    %saveas(gcf, sprintf('Fig_Pulse_%s_time_domain.png', shape_file));

    %% Plot 2: input FFT with accepted/rejected bins
    figure;
    semilogy(f_pos(idx_plot), I_mag(idx_plot), 'LineWidth', 1); hold on;

    acc_plot = idx_plot & accepted;
    rej_plot = idx_plot & rejected;

    semilogy(f_pos(acc_plot), I_mag(acc_plot), 'bo', 'MarkerSize', 4);
    semilogy(f_pos(rej_plot), I_mag(rej_plot), 'rx', 'MarkerSize', 4);
    yline(threshold, '--', 'Threshold = 1% max');

    grid on;
    xlabel('Frequency (Hz)');
    ylabel('|I(f)|');
    title(sprintf('%s pulse input FFT with accepted/rejected bins', strrep(shape_file, '_', '-')));
    legend('|I(f)|', 'Accepted bins', 'Rejected bins', 'Threshold', 'Location', 'best');
    %saveas(gcf, sprintf('Fig_Pulse_%s_accept_reject_fft.png', shape_file));

    %% Plot 3: impedance magnitude with accepted/rejected bins
    figure;
    plot(f_pos(idx_plot), abs(Z_true_pos(idx_plot)), 'k-', 'LineWidth', 1.5); hold on;
    plot(f_pos(acc_plot), abs(Z_est_pos(acc_plot)), 'bo', 'MarkerSize', 4);
    plot(f_pos(rej_plot), abs(Z_est_pos(rej_plot)), 'rx', 'MarkerSize', 4);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('|Z| (\Omega)');
    title(sprintf('%s pulse-estimated impedance magnitude', strrep(shape_file, '_', '-')));
    legend('Analytical R-L', 'Accepted estimate', 'Rejected estimate', 'Location', 'best');
    %saveas(gcf, sprintf('Fig_Pulse_%s_impedance_magnitude.png', shape_file));

    %% Plot 4: complex error with accepted/rejected bins
    figure;
    semilogy(f_pos(acc_plot), complex_error_all(acc_plot), 'bo', 'MarkerSize', 4); hold on;
    semilogy(f_pos(rej_plot), complex_error_all(rej_plot), 'rx', 'MarkerSize', 4);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Complex relative error (%)');
    title(sprintf('%s pulse impedance error with accepted/rejected bins', strrep(shape_file, '_', '-')));
    legend('Accepted bins', 'Rejected bins', 'Location', 'best');
    %saveas(gcf, sprintf('Fig_Pulse_%s_complex_error.png', shape_file));

    %% Save selected table to workspace with unique name
    assignin('base', sprintf('PulseTable_%s', shape_file), PulseTable);

end

%% Pulse-shape comparison summary
PulseShapeSummary = table(Shape, ...
    Accepted_Selected, Rejected_Selected, ...
    Accepted_0_to_50Hz, Rejected_0_to_50Hz, ...
    Mean_Complex_Error_Accepted_Selected, ...
    Best_Selected_Frequency_Hz, Best_Selected_Complex_Error_percent);

fprintf('\n=============================================\n');
fprintf('Pulse-shape comparison summary:\n');
disp(PulseShapeSummary);

%% Save summary table
writetable(PulseShapeSummary, 'PulseShapeSummary.csv');

%% Final recommendation based on this test
fprintf('\nSuggested interpretation:\n');
fprintf(['Compare the number of accepted selected points and the accepted-point error.\n' ...
    'Use the smoother pulse shape if it gives reliable low-frequency points and avoids severe weak-bin errors.\n' ...
    'Rejected bins should not be used for final impedance comparison.\n']);

%% Local function for phase wrapping
function angle_wrapped = wrapTo180_local(angle_deg)
    angle_wrapped = mod(angle_deg + 180, 360) - 180;
end