% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Performs the initial abc-waveform Richardson validation for the average GFMI model.

clear; clc; close all;
%% Repository paths
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

averageDataRoot = fullfile(repoRoot, 'data', 'student_a', 'average_gfmi');
op1Folder = fullfile(averageDataRoot, 'op1');
op2Folder = fullfile(averageDataRoot, 'op2');

resultRoot = fullfile(averageDataRoot, 'results');
figureFolder = fullfile(repoRoot, 'figures', 'student_a', 'average_gfmi');

if ~isfolder(op1Folder)
    error('OP1 data folder not found: %s', op1Folder);
end

if ~isfolder(op2Folder)
    error('OP2 data folder not found: %s', op2Folder);
end

if ~isfolder(resultRoot)
    mkdir(resultRoot);
end

if ~isfolder(figureFolder)
    mkdir(figureFolder);
end
%% Folder and files
dataFolder = op1Folder;

files.us50 = fullfile(dataFolder, 'GFMI_SCR2_50us.csv');
files.us25 = fullfile(dataFolder, 'GFMI_SCR2_25us.csv');
files.us12 = fullfile(dataFolder, 'GFMI_SCR2_12p5us.csv');
initialResultFolder = fullfile(resultRoot, 'initial_solverstep_validation');

if ~isfolder(initialResultFolder)
    mkdir(initialResultFolder);
end

%% Richardson settings
hCoarse = 50e-6;
hFine   = 25e-6;
p       = 2;
r       = hCoarse / hFine;
denom   = r^p - 1;

%% Analysis window
t_start_use = 2;
t_end_use   = 10;

%% Load data
D.us50 = load_csv(files.us50);
D.us25 = load_csv(files.us25);
D.us12 = load_csv(files.us12);

t = D.us50.Time;

%% Check time alignment
check_time_alignment(D.us50.Time, D.us25.Time, D.us12.Time);

%% Extract abc voltage and current
methods = {'us50','us25','us12'};

for k = 1:length(methods)
    m = methods{k};

    V.(m) = get_abc(D.(m), 'Vpoc');
    I.(m) = get_abc(D.(m), 'Ipoc');
end

%% Richardson correction in abc domain
V.RE = V.us25 + (V.us25 - V.us50) / denom;
I.RE = I.us25 + (I.us25 - I.us50) / denom;

%% Select steady-state window
idx = (t >= t_start_use) & (t < t_end_use);
t_use = t(idx);

%% Calculate waveform errors
[VoltageResults, V_rel] = calc_abc_errors(V, idx);
[CurrentResults, I_rel] = calc_abc_errors(I, idx);

fprintf('\nVoltage waveform error table compared with 12.5 us reference:\n');
disp(VoltageResults);

fprintf('\nCurrent waveform error table compared with 12.5 us reference:\n');
disp(CurrentResults);

%% Save results
writetable(VoltageResults, fullfile(initialResultFolder, 'GFMI_voltage_Richardson_results.csv'));
writetable(CurrentResults, fullfile(initialResultFolder, 'GFMI_current_Richardson_results.csv'));

%% Operating point check using 25 us case
P_mean    = mean(D.us25.Ppoc(idx));
Q_mean    = mean(D.us25.Qpoc(idx));
Vrms_mean = mean(D.us25.Vpoc_rms(idx));

fprintf('\nOperating point check using 25 us run:\n');
fprintf('Mean Ppoc     = %.6f\n', P_mean);
fprintf('Mean Qpoc     = %.6f\n', Q_mean);
fprintf('Mean Vpoc_rms = %.6f\n', Vrms_mean);

%% Plot relative RMS error bar charts
plot_error_bar(V_rel, 'GFMI PCC voltage relative RMS error compared with 12.5 us');
plot_error_bar(I_rel, 'GFMI PCC current relative RMS error compared with 12.5 us');

%% Waveform and error plots for phase A only
plot_phase_A(t_use, V, idx, 'V_{poc,A}', 'GFMI PCC voltage phase A');
plot_phase_A(t_use, I, idx, 'I_{poc,A}', 'GFMI PCC current phase A');

plot_phase_A_error(t_use, V, idx, 'Voltage error', 'GFMI voltage phase-A error compared with 12.5 us reference');
plot_phase_A_error(t_use, I, idx, 'Current error', 'GFMI current phase-A error compared with 12.5 us reference');

%% Local functions

function D = load_csv(file)
    D = readtable(file);
    D.Properties.VariableNames = strtrim(D.Properties.VariableNames);
end

function check_time_alignment(t50, t25, t12)
    if length(t50) ~= length(t25) || length(t50) ~= length(t12)
        error('Files do not have the same number of samples.');
    end
end

function Xabc = get_abc(D, prefix)
    Xabc = [D.([prefix '_A']), D.([prefix '_B']), D.([prefix '_C'])];
end

function [Results, RelErr] = calc_abc_errors(X, idx)
    phaseNames = {'A'; 'B'; 'C'};
    methodNames = {'us50','us25','RE'};
    methodLabels = ["50 us", "25 us", "Richardson"];
    MaxErr = zeros(3,3);
    RMSErr = zeros(3,3);
    RelErr = zeros(3,3);
    BestMethod = strings(3,1);
    Xref = X.us12(idx,:);
    for ph = 1:3
        ref_rms = rms(Xref(:,ph));
        for m = 1:3
            method = methodNames{m};
            err = X.(method)(idx,ph) - Xref(:,ph);
            MaxErr(ph,m) = max(abs(err));
            RMSErr(ph,m) = rms(err);
            RelErr(ph,m) = RMSErr(ph,m) / ref_rms * 100;
        end

        [~, bestIdx] = min(RelErr(ph,:));
        BestMethod(ph) = methodLabels(bestIdx);
    end

    Results = table( ...
        phaseNames, ...
        MaxErr(:,1), MaxErr(:,2), MaxErr(:,3), ...
        RMSErr(:,1), RMSErr(:,2), RMSErr(:,3), ...
        RelErr(:,1), RelErr(:,2), RelErr(:,3), ...
        BestMethod, ...
        'VariableNames', {'Phase', ...
        'MaxErr_50', 'MaxErr_25', 'MaxErr_RE', ...
        'RMSErr_50', 'RMSErr_25', 'RMSErr_RE', ...
        'RelRMSErr_50_percent', 'RelRMSErr_25_percent', 'RelRMSErr_RE_percent', ...
        'BestMethod'});
end

function plot_error_bar(RelErr, plotTitle)
    figure;
    bar(RelErr);
    grid on;
    xlabel('Phase');
    ylabel('Relative RMS error (%)');
    title(plotTitle);
    set(gca, 'XTickLabel', {'A','B','C'});
    legend('50 us', '25 us', 'Richardson', 'Location', 'best');
end

function plot_phase_A(t_use, X, idx, yLabelText, plotTitle)
    figure;
    plot(t_use, X.us50(idx,1), 'LineWidth', 1); hold on;
    plot(t_use, X.us25(idx,1), '--', 'LineWidth', 1);
    plot(t_use, X.RE(idx,1), ':', 'LineWidth', 1.8);
    plot(t_use, X.us12(idx,1), '-.', 'LineWidth', 1);
    grid on;
    xlabel('Time (s)');
    ylabel(yLabelText);
    title([plotTitle ': 50 us, 25 us, Richardson, 12.5 us']);
    legend('50 us', '25 us', 'Richardson', '12.5 us', 'Location', 'best');

end

function plot_phase_A_error(t_use, X, idx, yLabelText, plotTitle)
    figure;
    plot(t_use, abs(X.us50(idx,1) - X.us12(idx,1)), 'LineWidth', 1); hold on;
    plot(t_use, abs(X.us25(idx,1) - X.us12(idx,1)), '--', 'LineWidth', 1);
    plot(t_use, abs(X.RE(idx,1) - X.us12(idx,1)), ':', 'LineWidth', 1.8);
    grid on;
    xlabel('Time (s)');
    ylabel(yLabelText);
    title(plotTitle);
    legend('50 us error', '25 us error', 'Richardson error', 'Location', 'best');
end