% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares rectangular, raised-cosine and Gaussian pulses using analytical RL impedance error.
clear; clc; close all;

%% Fig. 2 - Pulse FFT Error Screening
% Compares rectangular, raised-cosine, and Gaussian pulse shapes
% using relative impedance error at selected low-frequency points.

%% Simulation settings
fs = 1000;                  % Sampling frequency [Hz]
Tsim = 20;                  % Simulation time [s]
t = 0:1/fs:Tsim-1/fs;       % Time vector
N = length(t);

f_fft = [0:N/2-1 -N/2:-1] * (fs/N);
w_fft = 2*pi*f_fft;

%% Analytical C0 R-L impedance
R = 1;                      % Ohm
L = 0.1;                    % H
Z_true = R + 1j*w_fft*L;

%% Pulse settings
pulse_amp = 1;
pulse_width = 0.2;
pulse_start = 2;
pulse_center = pulse_start + pulse_width/2;

pulse_names = {'Rectangular', 'Raised-cosine', 'Gaussian'};

%% Test settings
f_check = [0.2 0.5 1 2 5];      % Selected low-frequency points [Hz]
noise_level = 0.001;            % 0.1% voltage noise
accept_limit = 5;               % 5% relative error limit
weak_threshold_ratio = 0.01;    % 1% of maximum |I_fft|

rel_error_percent = zeros(length(pulse_names), length(f_check));
decision = strings(length(pulse_names), length(f_check));

%% Repeatable noise
rng(1);
noise_base = randn(size(t));

%% Run each pulse shape
for p = 1:length(pulse_names)

    I_time = zeros(size(t));

    switch pulse_names{p}

        case 'Rectangular'
            pulse_index = t >= pulse_start & t < pulse_start + pulse_width;
            I_time(pulse_index) = pulse_amp;

        case 'Raised-cosine'
            pulse_index = t >= pulse_start & t < pulse_start + pulse_width;
            tau = (t(pulse_index) - pulse_start) / pulse_width;
            I_time(pulse_index) = pulse_amp * 0.5 .* (1 - cos(2*pi*tau));

        case 'Gaussian'
            sigma = pulse_width/6;
            I_time = pulse_amp * exp(-((t - pulse_center).^2) / (2*sigma^2));
    end

    %% Analytical voltage response
    I_fft = fft(I_time);
    V_fft_clean = Z_true .* I_fft;
    V_time_clean = real(ifft(V_fft_clean));

    %% Add voltage noise
    V_time = V_time_clean + noise_level * max(abs(V_time_clean)) * noise_base;

    %% FFT impedance estimate
    V_fft_est = fft(V_time);
    Z_est = V_fft_est ./ I_fft;

    %% Positive frequencies only
    f_pos = f_fft(1:N/2);
    Z_true_pos = Z_true(1:N/2);
    Z_est_pos = Z_est(1:N/2);
    I_fft_pos = I_fft(1:N/2);

    %% Weak-current FFT screening
    threshold = weak_threshold_ratio * max(abs(I_fft_pos));
    valid_idx = abs(I_fft_pos) > threshold;

    %% Relative error at selected frequencies
    for k = 1:length(f_check)

        [~, idx_f] = min(abs(f_pos - f_check(k)));

        Zt = Z_true_pos(idx_f);
        Ze = Z_est_pos(idx_f);

        err = abs(Ze - Zt) / abs(Zt) * 100;
        rel_error_percent(p,k) = err;

        if ~valid_idx(idx_f)
            decision(p,k) = "Reject - weak FFT";
        elseif err < accept_limit
            decision(p,k) = "Accept";
        else
            decision(p,k) = "Reject - high error";
        end
    end
end

%% Print results
fprintf('\nPulse FFT error screening results:\n');
fprintf('Pulse width = %.2f s, Tsim = %.1f s, Noise level = %.2f%%\n\n', ...
    pulse_width, Tsim, noise_level*100);

for p = 1:length(pulse_names)
    fprintf('%s:\n', pulse_names{p});
    for k = 1:length(f_check)
        fprintf('  %.1f Hz: error = %.3g%%, %s\n', ...
            f_check(k), rel_error_percent(p,k), decision(p,k));
    end
    fprintf('\n');
end

%% Cap very large errors for plotting only
plot_cap = 10;
error_plot = rel_error_percent;
error_plot(error_plot > plot_cap) = plot_cap;

%% Final Fig. 2 plot
figure('Units','centimeters','Position',[2 2 9.5 7]);

semilogx(f_check, error_plot(1,:), '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 6); hold on;

semilogx(f_check, error_plot(2,:), '-s', ...
    'LineWidth', 1.5, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'w');

semilogx(f_check, error_plot(3,:), '-^', ...
    'LineWidth', 1.5, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'w');

yline(accept_limit, '--', '5% limit', ...
    'LineWidth', 1.2, ...
    'FontSize', 7, ...
    'LabelHorizontalAlignment', 'left');

grid on;
box on;

xlabel('Frequency (Hz)', 'FontSize', 9);
ylabel('Relative Error (%)', 'FontSize', 9);
title('Pulse FFT Error Screening', 'FontSize', 10, 'FontWeight', 'bold');

set(gca, 'XTick', f_check);
set(gca, 'XTickLabel', {'0.2','0.5','1','2','5'});

xlim([0.18 5.5]);
ylim([0 10.5]);

legend('Rectangular', 'Raised-cosine', 'Gaussian', '5% limit', ...
    'Location', 'northwest', ...
    'FontSize', 7, ...
    'Box', 'on');

set(gca, 'FontSize', 8);
set(gcf, 'Color', 'w');

set(gca, 'LooseInset', max(get(gca,'TightInset'), 0.05));

exportgraphics(gcf, 'fig2_pulse_error_screening_final.png', ...
    'Resolution', 600, ...
    'ContentType', 'image');