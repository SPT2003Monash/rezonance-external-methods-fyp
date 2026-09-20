% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Validates solver-timestep Richardson extrapolation for a PSCAD RL circuit at multiple frequencies.

clear; clc; close all;
%% Repository paths
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

dataFolder = fullfile(repoRoot, 'data', 'student_a', 'analytical_rlc');
figureFolder = fullfile(repoRoot, 'figures', 'student_a', 'analytical_rlc');
resultFolder = fullfile(repoRoot, 'data', 'student_a','analytical_rlc', 'results');

if ~isfolder(dataFolder)
    error('Analytical data folder not found: %s', dataFolder);
end

if ~isfolder(figureFolder)
    mkdir(figureFolder);
end

if ~isfolder(resultFolder)
    mkdir(resultFolder);
end

%% Circuit parameters
R = 1;          % ohm
L = 0.1;        % H

%% Frequencies tested
freqs = [1 5 10 20];   % Hz

%% Richardson settings
h0 = 50e-6;
h1 = 25e-6;
p  = 2;
r  = h0/h1;

%% Time window for phasor extraction
t_start_use = 2;
t_end_use   = 10;

%% Storage arrays
Z50_all   = zeros(length(freqs),1);
Z25_all   = zeros(length(freqs),1);
ZRE_all   = zeros(length(freqs),1);
Z12_all   = zeros(length(freqs),1);
Ztrue_all = zeros(length(freqs),1);

mag_err50 = zeros(length(freqs),1);
mag_err25 = zeros(length(freqs),1);
mag_errRE = zeros(length(freqs),1);
mag_err12 = zeros(length(freqs),1);

phase_err50 = zeros(length(freqs),1);
phase_err25 = zeros(length(freqs),1);
phase_errRE = zeros(length(freqs),1);
phase_err12 = zeros(length(freqs),1);

complex_err50 = zeros(length(freqs),1);
complex_err25 = zeros(length(freqs),1);
complex_errRE = zeros(length(freqs),1);
complex_err12 = zeros(length(freqs),1);

%% Loop over each frequency
for k = 1:length(freqs)

    f0 = freqs(k);

    fprintf('\n=============================================\n');
    fprintf('Processing source frequency = %g Hz\n', f0);

    %% File names
    file50 = sprintf('RL_f%g_50us.csv', f0);
    file25 = sprintf('RL_f%g_25us.csv', f0);
    file12 = sprintf('RL_f%g_12p5us.csv', f0);

    %% Load data
    data50 = readtable(fullfile(dataFolder, file50));
    data25 = readtable(fullfile(dataFolder, file25));
    data12 = readtable(fullfile(dataFolder, file12));
    

    data50.Properties.VariableNames = strtrim(data50.Properties.VariableNames);
    data25.Properties.VariableNames = strtrim(data25.Properties.VariableNames);
    data12.Properties.VariableNames = strtrim(data12.Properties.VariableNames);

    t50 = data50.Domain;
    t25 = data25.Domain;
    t12 = data12.Domain;

    i50 = data50.Ia;
    i25 = data25.Ia;
    i12 = data12.Ia;

    % Voltage sign correction
    v50 = -data50.Ea;
    v25 = -data25.Ea;
    v12 = -data12.Ea;

    %% Check time alignment
    if length(t50) ~= length(t25) || length(t50) ~= length(t12)
        error('Files for %g Hz do not have the same number of samples.', f0);
    end

    time_mismatch_25 = max(abs(t50 - t25));
    time_mismatch_12 = max(abs(t50 - t12));

    fprintf('Maximum time mismatch 50 vs 25 = %.3e s\n', time_mismatch_25);
    fprintf('Maximum time mismatch 50 vs 12.5 = %.3e s\n', time_mismatch_12);

    if time_mismatch_25 > 1e-12 || time_mismatch_12 > 1e-12
        error('Time vectors are not aligned for %g Hz.', f0);
    end

    t = t50;

    %% Richardson correction in time domain
    vRE = v25 + (v25 - v50)/(r^p - 1);
    iRE = i25 + (i25 - i50)/(r^p - 1);

    %% Select steady-state section
    idx = (t >= t_start_use) & (t <= t_end_use);

    t_use = t(idx);

    v50_use = v50(idx);
    v25_use = v25(idx);
    vRE_use = vRE(idx);
    v12_use = v12(idx);

    i50_use = i50(idx);
    i25_use = i25(idx);
    iRE_use = iRE(idx);
    i12_use = i12(idx);

    %% Time-domain current plot for 1 Hz case only
    if f0 == 1
        figure;
        plot(t_use, i50_use, 'LineWidth', 1); hold on;
        plot(t_use, i25_use, '--', 'LineWidth', 1);
        plot(t_use, iRE_use, ':', 'LineWidth', 1.8);
        plot(t_use, i12_use, '-.', 'LineWidth', 1);
        grid on;

        xlabel('Time (s)');
        ylabel('Current I_a (A)');
        title('Time-domain current: 50 us, 25 us, Richardson, and 12.5 us reference');
        legend('i_{50}', 'i_{25}', 'i_{RE}', 'i_{12.5}', 'Location', 'best');

        saveas(gcf, 'Fig_1_Time_domain_current_Richardson.png');
    end
    %% Estimate impedance
    Z50 = estimate_Z_phasor(t_use, v50_use, i50_use, f0);
    Z25 = estimate_Z_phasor(t_use, v25_use, i25_use, f0);
    ZRE = estimate_Z_phasor(t_use, vRE_use, iRE_use, f0);
    Z12 = estimate_Z_phasor(t_use, v12_use, i12_use, f0);

    %% Analytical true impedance
    Ztrue = R + 1j*2*pi*f0*L;

    %% Store impedance values
    Z50_all(k)   = Z50;
    Z25_all(k)   = Z25;
    ZRE_all(k)   = ZRE;
    Z12_all(k)   = Z12;
    Ztrue_all(k) = Ztrue;

    %% Calculate errors compared with analytical reference
    mag_true = abs(Ztrue);
    phase_true = angle(Ztrue)*180/pi;

    mag_err50(k) = abs(abs(Z50) - mag_true)/mag_true*100;
    mag_err25(k) = abs(abs(Z25) - mag_true)/mag_true*100;
    mag_errRE(k) = abs(abs(ZRE) - mag_true)/mag_true*100;
    mag_err12(k) = abs(abs(Z12) - mag_true)/mag_true*100;

    phase_err50(k) = abs(angle(Z50)*180/pi - phase_true);
    phase_err25(k) = abs(angle(Z25)*180/pi - phase_true);
    phase_errRE(k) = abs(angle(ZRE)*180/pi - phase_true);
    phase_err12(k) = abs(angle(Z12)*180/pi - phase_true);

    complex_err50(k) = abs(Z50 - Ztrue)/abs(Ztrue)*100;
    complex_err25(k) = abs(Z25 - Ztrue)/abs(Ztrue)*100;
    complex_errRE(k) = abs(ZRE - Ztrue)/abs(Ztrue)*100;
    complex_err12(k) = abs(Z12 - Ztrue)/abs(Ztrue)*100;

    %% Print result for this frequency
    fprintf('Ztrue = %.9f + j%.9f ohm\n', real(Ztrue), imag(Ztrue));
    fprintf('Z50   = %.9f + j%.9f ohm\n', real(Z50), imag(Z50));
    fprintf('Z25   = %.9f + j%.9f ohm\n', real(Z25), imag(Z25));
    fprintf('ZRE   = %.9f + j%.9f ohm\n', real(ZRE), imag(ZRE));
    fprintf('Z12.5 = %.9f + j%.9f ohm\n', real(Z12), imag(Z12));

end

%% Results table
ResultsTable = table(freqs(:), ...
    abs(Ztrue_all), abs(Z50_all), abs(Z25_all), abs(ZRE_all), abs(Z12_all), ...
    angle(Ztrue_all)*180/pi, angle(Z50_all)*180/pi, angle(Z25_all)*180/pi, ...
    angle(ZRE_all)*180/pi, angle(Z12_all)*180/pi, ...
    mag_err50, mag_err25, mag_errRE, mag_err12, ...
    phase_err50, phase_err25, phase_errRE, phase_err12, ...
    complex_err50, complex_err25, complex_errRE, complex_err12, ...
    'VariableNames', {'Frequency_Hz', ...
    'Ztrue_mag','Z50_mag','Z25_mag','ZRE_mag','Z12p5_mag', ...
    'Ztrue_phase_deg','Z50_phase_deg','Z25_phase_deg','ZRE_phase_deg','Z12p5_phase_deg', ...
    'MagErr_50_percent','MagErr_25_percent','MagErr_RE_percent','MagErr_12p5_percent', ...
    'PhaseErr_50_deg','PhaseErr_25_deg','PhaseErr_RE_deg','PhaseErr_12p5_deg', ...
    'ComplexErr_50_percent','ComplexErr_25_percent','ComplexErr_RE_percent','ComplexErr_12p5_percent'});

fprintf('\nFinal multi-frequency results table:\n');
disp(ResultsTable);

%% Plot: impedance magnitude
figure;
semilogx(freqs, abs(Ztrue_all), 'k-o', 'LineWidth', 1.5); hold on;
semilogx(freqs, abs(Z50_all), 'o-', 'LineWidth', 1);
semilogx(freqs, abs(Z25_all), 's--', 'LineWidth', 1);
semilogx(freqs, abs(ZRE_all), 'd:', 'LineWidth', 1.8);
semilogx(freqs, abs(Z12_all), 'x-.', 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('|Z| (\Omega)');
title('R-L impedance magnitude');
legend('Analytical', '50 us', '25 us', 'Richardson', '12.5 us', 'Location', 'best');

%% Plot: impedance phase
figure;
semilogx(freqs, angle(Ztrue_all)*180/pi, 'k-o', 'LineWidth', 1.5); hold on;
semilogx(freqs, angle(Z50_all)*180/pi, 'o-', 'LineWidth', 1);
semilogx(freqs, angle(Z25_all)*180/pi, 's--', 'LineWidth', 1);
semilogx(freqs, angle(ZRE_all)*180/pi, 'd:', 'LineWidth', 1.8);
semilogx(freqs, angle(Z12_all)*180/pi, 'x-.', 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('R-L impedance phase');
legend('Analytical', '50 us', '25 us', 'Richardson', '12.5 us', 'Location', 'best');

%% Plot: Bode magnitude error
figure;
semilogx(freqs, mag_err50, 'o-', 'LineWidth', 1); hold on;
semilogx(freqs, mag_err25, 's--', 'LineWidth', 1);
semilogx(freqs, mag_errRE, 'd:', 'LineWidth', 1.8);
semilogx(freqs, mag_err12, 'x-.', 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude error (%)');
title('Bode magnitude error compared with analytical R-L impedance');
legend('50 us', '25 us', 'Richardson', '12.5 us', 'Location', 'best');

%% Plot: Bode phase error
figure;
semilogx(freqs, phase_err50, 'o-', 'LineWidth', 1); hold on;
semilogx(freqs, phase_err25, 's--', 'LineWidth', 1);
semilogx(freqs, phase_errRE, 'd:', 'LineWidth', 1.8);
semilogx(freqs, phase_err12, 'x-.', 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase error (degrees)');
title('Bode phase error compared with analytical R-L impedance');
legend('50 us', '25 us', 'Richardson', '12.5 us', 'Location', 'best');

%% Plot: complex relative error
figure;
semilogx(freqs, complex_err50, 'o-', 'LineWidth', 1); hold on;
semilogx(freqs, complex_err25, 's--', 'LineWidth', 1);
semilogx(freqs, complex_errRE, 'd:', 'LineWidth', 1.8);
semilogx(freqs, complex_err12, 'x-.', 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('Complex relative error (%)');
title('Complex impedance error compared with analytical R-L impedance');
legend('50 us', '25 us', 'Richardson', '12.5 us', 'Location', 'best');

%% Local function
function Z = estimate_Z_phasor(t, v, i, f)

    v = v - mean(v);
    i = i - mean(i);

    basis = exp(-1j*2*pi*f*t);

    V = (2/length(t)) * sum(v .* basis);
    I = (2/length(t)) * sum(i .* basis);

    Z = V/I;
end