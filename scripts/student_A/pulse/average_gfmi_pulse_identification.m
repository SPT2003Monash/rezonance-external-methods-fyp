% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Performs the final average-GFMI pulse identification, screening, repeatability and linearity analysis.

clear; clc; close all;
%% Repository paths
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

dataFolder = fullfile(repoRoot, ...
    'data', 'student_a', 'pulse');

figureFolder = fullfile(repoRoot, ...
    'figures', 'student_a', 'pulse');

if ~isfolder(dataFolder)
    error('Pulse data folder not found:\n%s', dataFolder);
end

if ~isfolder(figureFolder)
    mkdir(figureFolder);
end
totalPulseTimer = tic;

%% USER SETTINGS


f_nom = 50;
w_nom = 2*pi*f_nom;

t_pre_start = 4.0;
t_pre_end = 4.9;
t_post_start = 4.0;
t_post_end = 24.0;
t_pulse_start = 2.0;
T_w = 0.5;

threshold_dB = -40;
cond_threshold = 10;
ref_tol_pct = 10;
screen_tol_pct = 20;
lin_tol_pct = 10;
min_cycles = 5;
f_min_plot = 0.1;
f_max_plot = 20;

f_ref_tones = [0.5, 1.0, 2.0, 3.0, 5.0];
ref_d_files = cellfun(@(x) fullfile(dataFolder,x), ...
    {'AVG_0p5Hz_dpert.csv','AVG_1Hz_dpert.csv','AVG_2Hz_dpert.csv', ...
     'AVG_3Hz_dpert.csv','AVG_5Hz_dpert.csv'}, 'UniformOutput', false);
ref_q_files = cellfun(@(x) fullfile(dataFolder,x), ...
    {'AVG_0p5Hz_qpert.csv','AVG_1Hz_qpert.csv','AVG_2Hz_qpert.csv', ...
     'AVG_3Hz_qpert.csv','AVG_5Hz_qpert.csv'}, 'UniformOutput', false);

%% EXTRACT FULL AND HALF AMPLITUDE PULSE ADMITTANCE

mainPulseTimer = tic;
fprintf('=== Full amplitude pulse ===\n');
[Ydd, Yqd, Ydq, Yqq, f_eval, accepted, cond_num] = pulse_admittance( ...
    fullfile(dataFolder,'AVG_dpulse_0p005pu.csv'), ...
    fullfile(dataFolder,'AVG_qpulse_0p5deg.csv'), ...
    f_nom, w_nom, t_pre_start, t_pre_end, t_post_start, t_post_end, ...
    threshold_dB, cond_threshold, min_cycles, f_min_plot, f_max_plot, true);
mainPulseMatlabTime = toc(mainPulseTimer);
fprintf('\nMAIN pulse MATLAB post-processing time = %.4f seconds\n', mainPulseMatlabTime);

%% REPEATABILITY CHECK: ORIGINAL FULL PULSE VS REPEATED FULL PULSE

fprintf('\n=== Full-pulse repeatability check ===\n');
[Ydd_r, Yqd_r, Ydq_r, Yqq_r, f_eval_r, accepted_r, ~] = pulse_admittance( ...
    fullfile(dataFolder,'AVG_dpulse_repeated.csv'), ...
    fullfile(dataFolder,'AVG_qpulse_repeated.csv'), ...
    f_nom, w_nom, t_pre_start, t_pre_end, t_post_start, t_post_end, ...
    threshold_dB, cond_threshold, min_cycles, f_min_plot, f_max_plot, false);

assert(length(f_eval_r) == length(f_eval), ...
    'Original and repeat frequency grids have different lengths.');
assert(max(abs(f_eval_r - f_eval)) < 1e-9, ...
    'Original and repeat frequency values do not match.');

common_bins = accepted & accepted_r;
union_bins = accepted | accepted_r;
bin_agreement_pct = 100 * sum(common_bins) / max(sum(union_bins),1);

fprintf('Original qualified bins: %d\n', sum(accepted));
fprintf('Repeat qualified bins:   %d\n', sum(accepted_r));
fprintf('Common qualified bins:   %d\n', sum(common_bins));
fprintf('Bin agreement:            %.2f%%\n', bin_agreement_pct);

repeat_names = {'Ydd','Yqd','Ydq','Yqq'};
Y_original = {Ydd, Yqd, Ydq, Yqq};
Y_repeat = {Ydd_r, Yqd_r, Ydq_r, Yqq_r};

fprintf('\nElement  MedianErr%%   P95Err%%    MaxErr%%    MaxFreq(Hz)\n');
fprintf('---------------------------------------------------------\n');
for k = 1:4
    yo = Y_original{k};
    yr = Y_repeat{k};
    valid = common_bins & isfinite(real(yo)) & isfinite(imag(yo)) & ...
            isfinite(real(yr)) & isfinite(imag(yr));
    err = 100 * abs(yr(valid) - yo(valid)) ./ max(abs(yo(valid)),1e-12);
    f_valid = f_eval(valid);
    median_err = median(err);
    p95_err = prctile(err,95);
    [max_err,max_idx] = max(err);
    fprintf('%-7s  %10.6f  %10.6f  %10.6f  %10.4f\n', ...
        repeat_names{k}, median_err, p95_err, max_err, f_valid(max_idx));
end

fprintf('\n=== Half amplitude pulse (linearity check) ===\n');
[Ydd_h, Yqd_h, Ydq_h, Yqq_h, f_eval_h, accepted_h, ~] = pulse_admittance( ...
    fullfile(dataFolder,'AVG_dpulse_0p0025pu.csv'), ...
    fullfile(dataFolder,'AVG_qpulse_0p25deg.csv'), ...
    f_nom, w_nom, t_pre_start, t_pre_end, t_post_start, t_post_end, ...
    threshold_dB, cond_threshold, min_cycles, f_min_plot, f_max_plot, true);

rejected = ~accepted;
elements = {'Ydd','Yqd','Ydq','Yqq'};
Y_full = {Ydd, Yqd, Ydq, Yqq};
Y_half = {Ydd_h, Yqd_h, Ydq_h, Yqq_h};
names_Y = {'Y_{dd}','Y_{qd}','Y_{dq}','Y_{qq}'};

%% SINUSOIDAL REFERENCE EXTRACTION

fprintf('\n=== Sinusoidal reference admittance ===\n');
Ydd_ref = NaN(length(f_ref_tones),1);
Yqd_ref = NaN(length(f_ref_tones),1);
Ydq_ref = NaN(length(f_ref_tones),1);
Yqq_ref = NaN(length(f_ref_tones),1);
have_ref = false;

for fi = 1:length(f_ref_tones)
    fr = f_ref_tones(fi);
    if ~isfile(ref_d_files{fi}) || ~isfile(ref_q_files{fi})
        warning('Missing reference files for %.1f Hz', fr);
        continue;
    end
    Ymat = singletone_admittance(ref_d_files{fi}, ref_q_files{fi}, fr, f_nom, t_pulse_start);
    Ydd_ref(fi) = Ymat(1,1);
    Yqd_ref(fi) = Ymat(2,1);
    Ydq_ref(fi) = Ymat(1,2);
    Yqq_ref(fi) = Ymat(2,2);
    fprintf('%.1f Hz: |Ydd|=%.4e |Yqd|=%.4e |Ydq|=%.4e |Yqq|=%.4e\n', ...
        fr, abs(Ydd_ref(fi)), abs(Yqd_ref(fi)), abs(Ydq_ref(fi)), abs(Yqq_ref(fi)));
    have_ref = true;
end
Y_ref = {Ydd_ref, Yqd_ref, Ydq_ref, Yqq_ref};

%% REFERENCE COMPARISON TABLE

if have_ref
    fprintf('\n--- Pulse vs Single-Tone Reference ---\n');
    fprintf('%-6s %-6s %-12s %-12s %-10s %-8s %-22s\n', ...
        'Freq','Elem','Pulse |Y|','Ref |Y|','Err%','Cycles','Decision');
    fprintf('%s\n', repmat('-',1,80));

    for fi = 1:length(f_ref_tones)
        fr = f_ref_tones(fi);
        if isnan(Ydd_ref(fi))
            continue;
        end
        [~,ki] = min(abs(f_eval - fr));
        n_cycles = fr * (t_post_end - t_post_start);
        Yp = [Ydd(ki), Yqd(ki), Ydq(ki), Yqq(ki)];
        Yr = [Ydd_ref(fi), Yqd_ref(fi), Ydq_ref(fi), Yqq_ref(fi)];

        for ei = 1:4
            if isnan(Yp(ei)) || ~accepted(ki)
                dec = 'REJECTED (bin)';
                err = NaN;
            else
                err = abs(Yp(ei)-Yr(ei))/abs(Yr(ei))*100;
                if n_cycles < min_cycles
                    dec = 'SCREENING (cycles)';
                elseif err < ref_tol_pct
                    dec = 'ACCEPTED';
                elseif err < screen_tol_pct
                    dec = 'SCREENING ONLY';
                else
                    dec = 'REJECTED (error)';
                end
            end
            fprintf('%-6.1f %-6s %-12.4e %-12.4e %-10s %-8.1f %-22s\n', ...
                fr, elements{ei}, abs(Yp(ei)), abs(Yr(ei)), sprintf('%.2f%%',err), n_cycles, dec);
        end
        fprintf('\n');
    end
end

%% LINEARITY CHECK TABLE

fprintf('\n--- Linearity Check: Full vs Half Amplitude ---\n');
fprintf('%-6s %-6s %-14s %-14s %-10s %-18s\n', ...
    'Freq','Elem','|Y| full','|Y| half','Change%','Linear?');
fprintf('%s\n', repmat('-',1,65));

for fi = 1:length(f_ref_tones)
    fr = f_ref_tones(fi);
    [~,kf] = min(abs(f_eval - fr));
    [~,kh] = min(abs(f_eval_h - fr));
    Yf = [Ydd(kf), Yqd(kf), Ydq(kf), Yqq(kf)];
    Yh = [Ydd_h(kh), Yqd_h(kh), Ydq_h(kh), Yqq_h(kh)];

    for ei = 1:4
        if isnan(Yf(ei)) || isnan(Yh(ei))
            fprintf('%-6.1f %-6s %-14s %-14s %-10s %-18s\n', ...
                fr, elements{ei}, 'N/A', 'N/A', 'N/A', 'N/A (rejected)');
        else
            chg = abs(Yf(ei)-Yh(ei))/abs(Yf(ei))*100;
            lin = 'LINEAR OK';
            if chg >= lin_tol_pct
                lin = 'NONLINEAR WARNING';
            end
            fprintf('%-6.1f %-6s %-14.4e %-14.4e %-10s %-18s\n', ...
                fr, elements{ei}, abs(Yf(ei)), abs(Yh(ei)), sprintf('%.2f%%',chg), lin);
        end
    end
    fprintf('\n');
end

%% ACCEPTED BIN SUMMARY

fprintf('\n--- Accepted Bin Summary ---\n');
fprintf('%-10s %-14s %-14s %-14s %-14s %-8s %-8s\n', ...
    'Freq(Hz)','|Ydd|','|Yqd|','|Ydq|','|Yqq|','cond#','Cycles');
for k = 1:length(f_eval)
    if accepted(k)
        fprintf('%-10.4f %-14.4e %-14.4e %-14.4e %-14.4e %-8.2f %-8.1f\n', ...
            f_eval(k), abs(Ydd(k)), abs(Yqd(k)), abs(Ydq(k)), abs(Yqq(k)), ...
            cond_num(k), f_eval(k)*(t_post_end-t_post_start));
    end
end

%% RUNTIME SAVING SUMMARY

n_acc = sum(accepted);
n_sine_runs = 2 * n_acc;
n_pulse_runs = 2;

fprintf('\n--- Runtime Saving ---\n');
fprintf('Accepted pulse bins:         %d\n', n_acc);
fprintf('Equivalent sinusoidal runs:  %d (d+q per frequency)\n', n_sine_runs);
fprintf('Pulse runs needed:           %d\n', n_pulse_runs);
fprintf('Run reduction:               %.1fx fewer runs\n', n_sine_runs/n_pulse_runs);

totalPulseMatlabTime_noPlots = toc(totalPulseTimer);
fprintf('\nTOTAL pulse MATLAB validation time excluding plots = %.4f seconds\n', totalPulseMatlabTime_noPlots);

%% PLOTS

% Figure 1: Time-domain (load signals from full amplitude files for plotting)
data_d_plot = readtable(fullfile(dataFolder,'AVG_dpulse_0p005pu.csv'));
data_q_plot = readtable(fullfile(dataFolder,'AVG_qpulse_0p5deg.csv'));
data_d_plot.Properties.VariableNames = strtrim(data_d_plot.Properties.VariableNames);
data_q_plot.Properties.VariableNames = strtrim(data_q_plot.Properties.VariableNames);
t_plot = data_d_plot.Time;

% Quick dq for plotting only
idx_pre_p = (t_plot >= t_pre_start) & (t_plot <= t_pre_end);

[Val_p, Vbe_p] = abc_to_alphabeta(data_q_plot.Vpoc_A, data_q_plot.Vpoc_B, data_q_plot.Vpoc_C);
th_p = angle(mean((Val_p(idx_pre_p)+1j*Vbe_p(idx_pre_p)) .* exp(-1j*w_nom*t_plot(idx_pre_p))));

[Valpha_dp, Vbeta_dp] = abc_to_alphabeta(data_d_plot.Vpoc_A, data_d_plot.Vpoc_B, data_d_plot.Vpoc_C);
sv_dp = Valpha_dp + 1j*Vbeta_dp;
th_dp = angle(mean(sv_dp(idx_pre_p) .* exp(-1j*w_nom*t_plot(idx_pre_p))));

theta_dp = w_nom*t_plot + th_dp;
theta_qp = w_nom*t_plot + th_p;

[Vd_dp,~] = abc_to_dq(data_d_plot.Vpoc_A, data_d_plot.Vpoc_B, data_d_plot.Vpoc_C, theta_dp);
[Id_dp,~] = abc_to_dq(data_d_plot.Ipoc_A, data_d_plot.Ipoc_B, data_d_plot.Ipoc_C, theta_dp);
[~,Vq_qp] = abc_to_dq(data_q_plot.Vpoc_A, data_q_plot.Vpoc_B, data_q_plot.Vpoc_C, theta_qp);
[~,Iq_qp] = abc_to_dq(data_q_plot.Ipoc_A, data_q_plot.Ipoc_B, data_q_plot.Ipoc_C, theta_qp);

dVd_dp = Vd_dp - mean(Vd_dp(idx_pre_p));
dId_dp = Id_dp - mean(Id_dp(idx_pre_p));
dVq_qp = Vq_qp - mean(Vq_qp(idx_pre_p));
dIq_qp = Iq_qp - mean(Iq_qp(idx_pre_p));

timePlotData = {dVd_dp, dId_dp, dVq_qp, dIq_qp};
timePlotColors = {'b','r','b','r'};
timePlotYLabels = {'\Delta V_d (pu)','\Delta I_d (pu)','\Delta V_q (pu)','\Delta I_q (pu)'};
timePlotTitles = {'d-pulse: \Delta V_d','d-pulse: \Delta I_d response', ...
                   'q-pulse: \Delta V_q','q-pulse: \Delta I_q response'};

figure('Name','Time-domain pulse signals');
for i = 1:4
    subplot(2,2,i)
    plot(t_plot, timePlotData{i}, timePlotColors{i}, 'LineWidth',1)
    xlabel('Time (s)')
    ylabel(timePlotYLabels{i})
    title(timePlotTitles{i})
    grid on
    xlim([1.5,8.5])
end

% Figures 2 and 3: Admittance magnitude and phase, each vs single-tone reference
plotAdmittanceSet('dq Admittance magnitude', f_eval, accepted, Y_full, names_Y, ...
    @(y) abs(y), '|%s| (S)', 'Admittance magnitude: ', have_ref, f_ref_tones, Y_ref);
plotAdmittanceSet('dq Admittance phase', f_eval, accepted, Y_full, names_Y, ...
    @(y) angle(y)*180/pi, 'Phase %s (deg)', 'Admittance phase: ', have_ref, f_ref_tones, Y_ref);

% Figure 4: Condition number
figure('Name','MIMO condition number');
semilogy(f_eval, cond_num, 'b', 'LineWidth',1.2)
hold on
semilogy(f_eval(accepted), cond_num(accepted), 'go', 'MarkerSize',5)
semilogy(f_eval(rejected), cond_num(rejected), 'rx', 'MarkerSize',5)
yline(cond_threshold, 'k--', sprintf('Threshold = %d', cond_threshold))
Twindow = t_post_end - t_post_start;
f_cycle_min = min_cycles/Twindow;
xline(f_cycle_min, 'm:', sprintf('%d-cycle limit = %.2f Hz', min_cycles, f_cycle_min), 'LineWidth', 1.2)
xlabel('Frequency (Hz)')
ylabel('Condition number \kappa(V_{dq})')
title('dq voltage matrix condition number')
legend('All bins','Accepted','Rejected','Threshold')
grid on

% Figure 5: Linearity
figure('Name','Linearity: Full vs Half Amplitude');
for ei = 1:4
    subplot(2,2,ei)
    plot(f_eval(accepted), abs(Y_full{ei}(accepted)), 'bo','MarkerSize',4)
    hold on
    plot(f_eval_h(accepted_h), abs(Y_half{ei}(accepted_h)), 'g^','MarkerSize',4)
    xlabel('Frequency (Hz)')
    ylabel(['|' names_Y{ei} '| (S)'])
    title(['Linearity: ' names_Y{ei}])
    legend('Full amplitude','Half amplitude','Location','best')
    grid on
end

% Figure 6: Error bar chart
if have_ref
    figure('Name','Pulse vs Reference Error');
    err_mat = zeros(4, length(f_ref_tones));
    for fi = 1:length(f_ref_tones)
        [~,ki] = min(abs(f_eval - f_ref_tones(fi)));
        for ei = 1:4
            Yp = Y_full{ei}(ki);
            Yr = Y_ref{ei}(fi);
            if isnan(Yp)
                err_mat(ei,fi) = NaN;
            else
                err_mat(ei,fi) = abs(Yp-Yr)/abs(Yr)*100;
            end
        end
    end
    bar(err_mat')
    set(gca,'XTickLabel',{'0.5 Hz','1 Hz','2 Hz','3 Hz'})
    ylabel('Error vs sinusoidal reference (%)')
    title('Pulse identification error at reference frequencies')
    legend('Ydd','Yqd','Ydq','Yqq','Location','best')
    yline(ref_tol_pct, 'k--', 'Accepted (10%)')
    yline(screen_tol_pct, 'r--', 'Screening (20%)')
    grid on
end

% Figure 7: Runtime saving
figure('Name','Runtime Saving');
bar([n_sine_runs, n_pulse_runs], 'FaceColor',[0.2 0.6 0.8])
set(gca,'XTickLabel',{sprintf('Sinusoidal\n(%d runs)',n_sine_runs), sprintf('Pulse\n(%d runs)',n_pulse_runs)})
ylabel('Number of PSCAD runs')
title(sprintf('Runs required for %d frequency points (%.1fx saving)', n_acc, n_sine_runs/n_pulse_runs))
grid on

%% SAVE RESULTS

f_acc = f_eval(accepted);
cycles_acc = f_acc * (t_post_end - t_post_start);
ResultsTable = table(f_acc, cycles_acc, ...
    real(Ydd(accepted)), imag(Ydd(accepted)), abs(Ydd(accepted)), angle(Ydd(accepted))*180/pi, ...
    real(Yqd(accepted)), imag(Yqd(accepted)), abs(Yqd(accepted)), angle(Yqd(accepted))*180/pi, ...
    real(Ydq(accepted)), imag(Ydq(accepted)), abs(Ydq(accepted)), angle(Ydq(accepted))*180/pi, ...
    real(Yqq(accepted)), imag(Yqq(accepted)), abs(Yqq(accepted)), angle(Yqq(accepted))*180/pi, ...
    cond_num(accepted), ...
    'VariableNames', {'Freq_Hz','Cycles_in_Window', ...
    'Re_Ydd','Im_Ydd','Mag_Ydd','Phase_Ydd_deg', ...
    'Re_Yqd','Im_Yqd','Mag_Yqd','Phase_Yqd_deg', ...
    'Re_Ydq','Im_Ydq','Mag_Ydq','Phase_Ydq_deg', ...
    'Re_Yqq','Im_Yqq','Mag_Yqq','Phase_Yqq_deg','CondNum'});
writetable(ResultsTable, fullfile(dataFolder,'AVG_pulse_admittance_accepted_24s.csv'));
fprintf('\nResults saved.\n');

%% ---------------------------------------------------------------
%  FUNCTIONS
%% ---------------------------------------------------------------

function plotAdmittanceSet(figName, f_eval, accepted, Yset, names, transformFcn, unitLabel, titlePrefix, haveRef, fRef, Yref)
    % Shared plotting routine for the magnitude and phase admittance figures.
    figure('Name', figName);
    for ei = 1:4
        subplot(2,2,ei)
        plot(f_eval(accepted), transformFcn(Yset{ei}(accepted)), 'bo', 'MarkerSize',5)
        hold on
        if haveRef
            plot(fRef, transformFcn(Yref{ei}), 'rs', 'MarkerSize',8, 'LineWidth',1.5)
            legend('Pulse (accepted)','Single-tone ref','Location','best')
        end
        xlabel('Frequency (Hz)')
        ylabel(sprintf(unitLabel, names{ei}))
        title([titlePrefix names{ei}])
        grid on
    end
end

function [Ydd,Yqd,Ydq,Yqq,f_eval,accepted,cond_num] = pulse_admittance(...
    file_d, file_q, f_nom, w_nom, ...
    t_pre_start, t_pre_end, t_post_start, t_post_end, ...
    threshold_dB, cond_threshold, min_cycles, ...
    f_min, f_max, verbose)

    dd = readtable(file_d);
    dd.Properties.VariableNames = strtrim(dd.Properties.VariableNames);
    dq = readtable(file_q);
    dq.Properties.VariableNames = strtrim(dq.Properties.VariableNames);

    t = dd.Time;
    dt = mean(diff(t));
    fs = 1/dt;
    if verbose
        fprintf('dt=%.2f us, fs=%.0f Hz\n', dt*1e6, fs);
    end
    Vabc_d = [dd.Vpoc_A, dd.Vpoc_B, dd.Vpoc_C];
    Iabc_d = [dd.Ipoc_A, dd.Ipoc_B, dd.Ipoc_C];
    Vabc_q = [dq.Vpoc_A, dq.Vpoc_B, dq.Vpoc_C];
    Iabc_q = [dq.Ipoc_A, dq.Ipoc_B, dq.Ipoc_C];

    idx_pre = (t>=t_pre_start) & (t<=t_pre_end);

    [Va,Vb] = abc_to_alphabeta(Vabc_d(:,1),Vabc_d(:,2),Vabc_d(:,3));
    th_d = angle(mean((Va(idx_pre)+1j*Vb(idx_pre)).*exp(-1j*w_nom*t(idx_pre))));
    [Va,Vb] = abc_to_alphabeta(Vabc_q(:,1),Vabc_q(:,2),Vabc_q(:,3));
    th_q = angle(mean((Va(idx_pre)+1j*Vb(idx_pre)).*exp(-1j*w_nom*t(idx_pre))));
    if verbose
        fprintf('d-offset=%.6f rad, q-offset=%.6f rad, diff=%.6f rad\n', th_d, th_q, abs(th_d-th_q));
    end

    theta_d = w_nom*t+th_d;
    theta_q = w_nom*t+th_q;

    [Vdd,Vqd] = abc_to_dq(Vabc_d(:,1),Vabc_d(:,2),Vabc_d(:,3),theta_d);
    [Idd,Iqd] = abc_to_dq(Iabc_d(:,1),Iabc_d(:,2),Iabc_d(:,3),theta_d);
    [Vdq,Vqq] = abc_to_dq(Vabc_q(:,1),Vabc_q(:,2),Vabc_q(:,3),theta_q);
    [Idq,Iqq] = abc_to_dq(Iabc_q(:,1),Iabc_q(:,2),Iabc_q(:,3),theta_q);
    dVd_d = Vdd-mean(Vdd(idx_pre));
    dVq_d = Vqd-mean(Vqd(idx_pre));
    dId_d = Idd-mean(Idd(idx_pre));
    dIq_d = Iqd-mean(Iqd(idx_pre));
    dVd_q = Vdq-mean(Vdq(idx_pre));
    dVq_q = Vqq-mean(Vqq(idx_pre));
    dId_q = Idq-mean(Idq(idx_pre));
    dIq_q = Iqq-mean(Iqq(idx_pre));

    idx_resp = (t>=t_post_start) & (t<t_post_end);
    N = sum(idx_resp);

    if verbose
        fprintf('Window: %.1f-%.1f s, N=%d\n', t_post_start, t_post_end, N);
        pre_amp = max(abs(dId_d(idx_pre)));
        late = max(abs(dId_d(end-round(0.1*N):end)));
        fprintf('Pre-pulse amp=%.3e, Late=%.3e\n', pre_amp, late);
        if late > 1.5*pre_amp
            warning('Response not decayed.');
        else
            fprintf('Decay check OK\n');
        end
    end

    % Rectangular FFT window.
    % The raised-cosine excitation is already smoothly tapered, and the
    % response begins and ends near the same settled operating point.
    fft_sig = @(x) fft(x(idx_resp))*2/N;

    F_Vd_d = fft_sig(dVd_d);
    F_Vq_d = fft_sig(dVq_d);
    F_Id_d = fft_sig(dId_d);
    F_Iq_d = fft_sig(dIq_d);
    F_Vd_q = fft_sig(dVd_q);
    F_Vq_q = fft_sig(dVq_q);
    F_Id_q = fft_sig(dId_q);
    F_Iq_q = fft_sig(dIq_q);
    f = (0:N-1)'*(fs/N);
    idx_p = f>=f_min & f<=f_max;
    f_eval = f(idx_p);
    sel = @(x) x(idx_p);
    F_Vd_d = sel(F_Vd_d);
    F_Vq_d = sel(F_Vq_d);
    F_Id_d = sel(F_Id_d);
    F_Iq_d = sel(F_Iq_d);
    F_Vd_q = sel(F_Vd_q);
    F_Vq_q = sel(F_Vq_q);
    F_Id_q = sel(F_Id_q);
    F_Iq_q = sel(F_Iq_q);
    % Minimum-cycle qualification
    Twindow = t_post_end - t_post_start;
    n_cycles = f_eval * Twindow;
    cycle_ok = n_cycles >= min_cycles;
    thr = 10^(threshold_dB/20);
    energy_ok = abs(F_Vd_d)>=thr*max(abs(F_Vd_d)) & abs(F_Vq_q)>=thr*max(abs(F_Vq_q));

    nf = length(f_eval);
    cond_num = zeros(nf,1);
    Yf = zeros(2,2,nf);
    for k = 1:nf
        V = [F_Vd_d(k),F_Vd_q(k); F_Vq_d(k),F_Vq_q(k)];
        I = [F_Id_d(k),F_Id_q(k); F_Iq_d(k),F_Iq_q(k)];
        cond_num(k) = cond(V);
        if cond_num(k) < cond_threshold*10
            Yf(:,:,k) = I/V;
        else
            Yf(:,:,k) = NaN(2,2);
        end
    end

    cond_ok = cond_num < cond_threshold;
    accepted = energy_ok & cond_ok & cycle_ok;
    Ydd = squeeze(Yf(1,1,:));
    Ydd(~accepted) = NaN;
    Yqd = squeeze(Yf(2,1,:));
    Yqd(~accepted) = NaN;
    Ydq = squeeze(Yf(1,2,:));
    Ydq(~accepted) = NaN;
    Yqq = squeeze(Yf(2,2,:));
    Yqq(~accepted) = NaN;

    if verbose
        fprintf('\n--- Broadband Bin Qualification ---\n');
        fprintf('FFT window duration: %.2f s\n', Twindow);
        fprintf('Minimum required cycles: %d\n', min_cycles);
        fprintf('Minimum cycle-qualified frequency: %.2f Hz\n', min_cycles/Twindow);
        fprintf('Energy-qualified bins: %d/%d\n', sum(energy_ok), nf);
        fprintf('Condition-qualified bins: %d/%d\n', sum(cond_ok), nf);
        fprintf('Cycle-qualified bins: %d/%d\n', sum(cycle_ok), nf);
        fprintf('Final accepted bins: %d/%d (%.1f%%)\n', sum(accepted), nf, 100*sum(accepted)/nf);
        fprintf('Bins rejected for insufficient cycles: %d\n', sum(~cycle_ok));
        if any(accepted)
            fprintf('Final accepted frequency range: %.2f to %.2f Hz\n', min(f_eval(accepted)), max(f_eval(accepted)));
        else
            warning('No frequency bins passed all qualification criteria.');
        end
    end
end

function Ymat = singletone_admittance(file_d, file_q, fr, f_nom, t_start)
    dd = readtable(file_d);
    dd.Properties.VariableNames = strtrim(dd.Properties.VariableNames);
    dq = readtable(file_q);
    dq.Properties.VariableNames = strtrim(dq.Properties.VariableNames);
    t = dd.Time;
    w_nom = 2*pi*f_nom;

    t1 = t_start+max(1.0,1/fr);
    nC = floor((max(t)-0.1-t1)*fr);
    if nC < 1
        idx = (t>=t_start+0.5) & (t<=max(t)-0.1);
    else
        idx = (t>=t1) & (t<t1+nC/fr);
    end
    t_use = t(idx);

    [Vdd,Vqd,Idd,Iqd] = to_dq_ss(dd,t,idx,w_nom);
    [Vdq,Vqq,Idq,Iqq] = to_dq_ss(dq,t,idx,w_nom);

    ph = @(x,f) estimate_phasor(t_use,x(idx),f);
    Vmat = [ph(Vdd,fr),ph(Vdq,fr); ph(Vqd,fr),ph(Vqq,fr)];
    Imat = [ph(Idd,fr),ph(Idq,fr); ph(Iqd,fr),ph(Iqq,fr)];
    Ymat = Imat/Vmat;
end

function [Vd,Vq,Id,Iq] = to_dq_ss(data, t, idx, w_nom)
    [Va,Vb] = abc_to_alphabeta(data.Vpoc_A,data.Vpoc_B,data.Vpoc_C);
    sv = Va+1j*Vb;
    th = angle(mean(sv(idx).*exp(-1j*w_nom*t(idx))));
    theta = w_nom*t+th;
    [Vd,Vq] = abc_to_dq(data.Vpoc_A,data.Vpoc_B,data.Vpoc_C,theta);
    [Id,Iq] = abc_to_dq(data.Ipoc_A,data.Ipoc_B,data.Ipoc_C,theta);
end

function [alpha,beta] = abc_to_alphabeta(a,b,c)
    alpha = sqrt(2/3)*(a-0.5*b-0.5*c);
    beta = sqrt(2/3)*((sqrt(3)/2)*b-(sqrt(3)/2)*c);
end

function [d,q] = abc_to_dq(a,b,c,theta)
    [alpha,beta] = abc_to_alphabeta(a,b,c);
    d = alpha.*cos(theta)+beta.*sin(theta);
    q = -alpha.*sin(theta)+beta.*cos(theta);
end

function X = estimate_phasor(t,x,f)
    x = x-mean(x);
    X = (2/length(t))*sum(x.*exp(-1j*2*pi*f*t));
end