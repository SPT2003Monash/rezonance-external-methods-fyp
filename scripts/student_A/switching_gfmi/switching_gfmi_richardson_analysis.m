% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Evaluates Richardson methods on envelope-extracted switching-GFMI waveforms.

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

dataFolder = fullfile(repoRoot, 'data', 'student_a', 'switching_gfmi');
resultFolder = fullfile(dataFolder, 'results');
figureFolder = fullfile(repoRoot, 'figures', 'student_a', 'switching_gfmi');

if ~isfolder(dataFolder)
    error('Switching-GFMI data folder not found: %s', dataFolder);
end
if ~isfolder(resultFolder); mkdir(resultFolder); end
if ~isfolder(figureFolder); mkdir(figureFolder); end
%% 1. SETTINGS ---------------------------------------------------------------

folder = dataFolder;
fpert = 10;       % perturbation frequency: 1, 5, or 10 Hz
f0    = 50;      % grid frequency [Hz]
fpwm  = 5000;    % PWM frequency [Hz]
fc    = 100;     % envelope LPF cutoff [Hz]

thetaWindow  = [0.5 1.9];   % clean pre-perturbation interval
phasorWindow = [3.0 10.0];  % admittance extraction interval

axesList = {'d','q'};
levels   = {'coarse','fine','finer','ref'};
stepTag = struct('coarse','12p5us','fine','6p25us','finer','3p125us','ref','0p1us');

fTag = sprintf('%dHz',fpert);
for a = 1:2
    for k = 1:numel(levels)
        ax = axesList{a}; lv = levels{k};
        files.(ax).(lv) = fullfile(folder, sprintf('Switching_GFMI_%s_%spert_%s.csv', fTag, ax, stepTag.(lv)));
    end
end

fprintf('\n=== SWITCHING GFMI — %.0f Hz FINAL RICHARDSON STUDY ===\n', fpert);

%% 2. LOAD + TIME-ALIGNMENT CHECK -------------------------------------------

for a = 1:2
    for k = 1:numel(levels)
        ax = axesList{a}; lv = levels{k};
        data.(ax).(lv) = loadCSV(files.(ax).(lv));
    end
end

t  = data.d.coarse.t;
dt = median(diff(t));
fs = 1/dt;
fprintf('Logging step = %.3f us | Sampling rate = %.0f Hz\n', dt*1e6, fs);

for a = 1:2
    for k = 1:numel(levels)
        ax = axesList{a}; lv = levels{k};
        tk = data.(ax).(lv).t;
        if length(tk) ~= length(t) || max(abs(tk-t)) > 1e-9
            error('%s.%s is not time-aligned.', ax, lv);
        end
    end
end
fprintf('All %d files are time aligned.\n', 2*numel(levels));

%% 3. PWM ENVELOPE EXTRACTION ------------------------------------------------

Nma = round((1/fpwm)/dt);
[bLP,aLP] = butter(4, fc/(fs/2), 'low');
fprintf('Moving average = %d samples (%.1f us) | LPF cutoff = %.1f Hz\n', Nma, Nma*dt*1e6, fc);

for a = 1:2
    for k = 1:numel(levels)
        ax = axesList{a}; lv = levels{k};
        data.(ax).(lv).Venv = getEnvelope(data.(ax).(lv).Vabc, Nma, bLP, aLP);
        data.(ax).(lv).Ienv = getEnvelope(data.(ax).(lv).Iabc, Nma, bLP, aLP);
    end
end

%% 4. COMMON dq REFERENCE ANGLE ----------------------------------------------

idxTheta = t >= thetaWindow(1) & t <= thetaWindow(2);
offD = thetaOffset(data.d.fine.Venv, t, idxTheta, 2*pi*f0);
offQ = thetaOffset(data.q.fine.Venv, t, idxTheta, 2*pi*f0);
offCommon = angle(exp(1j*offD) + exp(1j*offQ));
theta = 2*pi*f0*t + offCommon;

fprintf('\n=== COMMON dq ANGLE ===\n');
fprintf('6.25us d offset = %.9f rad (%.6f deg)\n', offD, rad2deg(offD));
fprintf('6.25us q offset = %.9f rad (%.6f deg)\n', offQ, rad2deg(offQ));
fprintf('Common offset   = %.9f rad (%.6f deg)\n', offCommon, rad2deg(offCommon));

%% 5. abc -> dq + RICHARDSON EXTRAPOLATION ------------------------------------
% p=1:   RE1   = 2*x6.25 - x12.5                          (cancels O(h))
% p=2:   RE2   = x6.25 + (x6.25-x12.5)/3                  (assumes p=2)
% mixed: mixed = (1/3)*x12.5 - 2*x6.25 + (8/3)*x3.125      (cancels O(h) and O(h^2))

signals = {'Vd','Vq','Id','Iq'};
for a = 1:2
    ax = axesList{a};
    for k = 1:numel(levels)
        lv = levels{k};
        S = data.(ax).(lv);
        [S.Vd,S.Vq,S.Id,S.Iq] = abc2dq(S.Venv, S.Ienv, theta);
        data.(ax).(lv) = S;
    end
    for s = 1:numel(signals)
        sig = signals{s};
        x12 = data.(ax).coarse.(sig);
        x6  = data.(ax).fine.(sig);
        x3  = data.(ax).finer.(sig);
        data.(ax).RE1.(sig)   = 2*x6 - x12;
        data.(ax).RE2.(sig)   = x6 + (x6-x12)/3;
        data.(ax).mixed.(sig) = (1/3)*x12 - 2*x6 + (8/3)*x3;
    end
end

%% 6. WAVEFORM ERRORS vs 0.1 us REFERENCE ------------------------------------

idxFit = t >= phasorWindow(1) & t < phasorWindow(2);
methods = {'coarse','fine','finer','RE1','RE2','mixed'};
labels  = {'12.5 us','6.25 us','3.125 us','p=1 RE','p=2 RE','Mixed RE'};

fprintf('\n=== dq ENVELOPE WAVEFORM ERRORS vs 0.1 us ===\n');
for a = 1:2
    ax = axesList{a};
    Vref = [data.(ax).ref.Vd(idxFit), data.(ax).ref.Vq(idxFit)];
    Iref = [data.(ax).ref.Id(idxFit), data.(ax).ref.Iq(idxFit)];

    fprintf('\n%s perturbation\n', upper(ax));
    for k = 1:numel(methods)
        lv = methods{k};
        Vnow = [data.(ax).(lv).Vd(idxFit), data.(ax).(lv).Vq(idxFit)];
        Inow = [data.(ax).(lv).Id(idxFit), data.(ax).(lv).Iq(idxFit)];
        fprintf('%-10s V = %.6f%%   I = %.6f%%\n', labels{k}, errMat(Vnow,Vref), errMat(Inow,Iref));
    end
end

%% 7. PHASORS + FULL 2x2 MATRICES --------------------------------------------

phasorLevels = {'coarse','fine','finer','RE1','RE2','mixed','ref'};
tFit = t(idxFit);

for a = 1:2
    ax = axesList{a};
    for k = 1:numel(phasorLevels)
        lv = phasorLevels{k};
        for s = 1:numel(signals)
            sig = signals{s};
            x = data.(ax).(lv).(sig);
            P.(ax).(lv).(sig) = phasorLS(tFit, x(idxFit), fpert);
        end
    end
end

fprintf('\n=== VOLTAGE MATRIX CONDITIONING ===\n');
for k = 1:numel(phasorLevels)
    lv = phasorLevels{k};
    Vmat.(lv) = [P.d.(lv).Vd, P.q.(lv).Vd; P.d.(lv).Vq, P.q.(lv).Vq];
    Imat.(lv) = [P.d.(lv).Id, P.q.(lv).Id; P.d.(lv).Iq, P.q.(lv).Iq];
    Yadm.(lv) = Imat.(lv) / Vmat.(lv);   % right division -- avoids inv(V)
    fprintf('%-7s cond(V) = %.6f\n', lv, cond(Vmat.(lv)));
end

%% 8. FULL-MATRIX AND ELEMENTWISE ERRORS --------------------------------------

matrixErr = zeros(1,numel(methods));
fprintf('\n=== COMPLETE dq ADMITTANCE MATRIX ERRORS ===\n');
for k = 1:numel(methods)
    matrixErr(k) = errMat(Yadm.(methods{k}), Yadm.ref);
    fprintf('%-10s = %.6f%%\n', labels{k}, matrixErr(k));
end
[bestErr,bestIdx] = min(matrixErr);
fprintf('\nBEST RESULT = %s (%.6f%%)\n', labels{bestIdx}, bestErr);

names = {'Ydd','Ydq','Yqd','Yqq'};
pos = [1 1; 1 2; 2 1; 2 2];
elemErr = zeros(4,numel(methods));

fprintf('\n=== ELEMENT-WISE COMPLEX ERRORS ===\n');
for e = 1:4
    r = pos(e,1); c = pos(e,2);
    fprintf('\n%s\n', names{e});
    for k = 1:numel(methods)
        elemErr(e,k) = errComplex(Yadm.(methods{k})(r,c), Yadm.ref(r,c));
        fprintf('%-10s = %.6f%%\n', labels{k}, elemErr(e,k));
    end
end

%% 9. OBSERVED CONVERGENCE (coarse->fine and fine->finer) ---------------------

pObs = zeros(4,1); ratioAngle = zeros(4,1);
pObs2 = zeros(4,1); ratioAngle2 = zeros(4,1);

fprintf('\n=== OBSERVED CONVERGENCE ===\n');
for e = 1:4
    r = pos(e,1); c = pos(e,2);
    ratio1 = (Yadm.coarse(r,c)-Yadm.ref(r,c)) / (Yadm.fine(r,c)-Yadm.ref(r,c));
    ratio2 = (Yadm.fine(r,c)-Yadm.ref(r,c))   / (Yadm.finer(r,c)-Yadm.ref(r,c));
    pObs(e) = log(abs(ratio1))/log(2);         ratioAngle(e)  = rad2deg(angle(ratio1));
    pObs2(e) = log(abs(ratio2))/log(2);        ratioAngle2(e) = rad2deg(angle(ratio2));
    fprintf('%s: 12.5->6.25 p=%.4f angle=%.4fdeg | 6.25->3.125 p=%.4f angle=%.4fdeg\n', ...
        names{e}, pObs(e), ratioAngle(e), pObs2(e), ratioAngle2(e));
end

%% 10. REFERENCE + SELECTED MATRIX --------------------------------------------

fprintf('\n=== 0.1 us REFERENCE MATRIX ===\n'); printYmatrix(Yadm.ref, names, pos);
bestMethod = methods{bestIdx};
fprintf('\n=== BEST METHOD MATRIX (%s) ===\n', labels{bestIdx}); printYmatrix(Yadm.(bestMethod), names, pos);

%% 11. SAVE RESULTS -------------------------------------------------------------

ResultTable = table(string(names)', elemErr(:,1),elemErr(:,2),elemErr(:,3),elemErr(:,4),elemErr(:,5),elemErr(:,6), ...
    pObs, ratioAngle, pObs2, ratioAngle2, ...
    'VariableNames',{'Element','Error_12p5us_pct','Error_6p25us_pct','Error_3p125us_pct', ...
    'Error_RE1_pct','Error_RE2_pct','Error_Mixed_pct', ...
    'ObservedP_12p5_to_6p25','AngleDeg_12p5_to_6p25','ObservedP_6p25_to_3p125','AngleDeg_6p25_to_3p125'});
disp(ResultTable);

outFile = fullfile(resultFolder,sprintf('SW_GFMI_%dHz_FINAL_Richardson_results.csv', fpert));writetable(ResultTable, outFile);
fprintf('\nSaved: %s\n', outFile);

%% 12. PLOTS --------------------------------------------------------------------

figure('Name','Full-matrix comparison');
bar(matrixErr); set(gca,'XTickLabel',labels);
ylabel('Full-matrix error (%)'); title(sprintf('%g Hz Switching GFMI — Richardson Comparison', fpert));
grid on; box on;

figure('Name','Element-wise comparison');
bar(elemErr); set(gca,'XTickLabel',names);
ylabel('Complex error (%)'); title(sprintf('%g Hz Switching GFMI — Element-wise Error', fpert));
legend(labels,'Location','best'); grid on; box on;

%% ============================ LOCAL FUNCTIONS ================================

function S = loadCSV(file)
    if ~isfile(file); error('File not found: %s', file); end
    T = readtable(file);
    T.Properties.VariableNames = strtrim(T.Properties.VariableNames);

    req = {'Time','Vpoc_A','Vpoc_B','Vpoc_C','Ipoc_A','Ipoc_B','Ipoc_C'};
    for k = 1:numel(req)
        if ~ismember(req{k}, T.Properties.VariableNames)
            error('Missing %s in %s', req{k}, file);
        end
    end
    S.t = T.Time;
    S.Vabc = [T.Vpoc_A, T.Vpoc_B, T.Vpoc_C];
    S.Iabc = [T.Ipoc_A, T.Ipoc_B, T.Ipoc_C];
end

function Xout = getEnvelope(Xin, N, b, a)
    Xout = zeros(size(Xin));
    for k = 1:3
        xMA = movmean(Xin(:,k), N, 'Endpoints','shrink');   % removes switching ripple
        Xout(:,k) = filtfilt(b, a, xMA);                     % zero-phase smoothing
    end
end

function off = thetaOffset(V, t, idx, w)
    [alpha,beta] = clarke(V(:,1), V(:,2), V(:,3));
    z = alpha + 1j*beta;
    off = angle(mean(z(idx) .* exp(-1j*w*t(idx))));
end

function [Vd,Vq,Id,Iq] = abc2dq(V, I, theta)
    [Va,Vb] = clarke(V(:,1), V(:,2), V(:,3));
    [Ia,Ib] = clarke(I(:,1), I(:,2), I(:,3));
    ct = cos(theta); st = sin(theta);
    Vd =  Va.*ct + Vb.*st;  Vq = -Va.*st + Vb.*ct;
    Id =  Ia.*ct + Ib.*st;  Iq = -Ia.*st + Ib.*ct;
end

function [alpha,beta] = clarke(a, b, c)
    alpha = sqrt(2/3)*(a - 0.5*b - 0.5*c);
    beta  = sqrt(2/3)*(sqrt(3)/2).*(b - c);
end

function X = phasorLS(t, x, f)
    t = t(:); x = x(:);
    tau = t - mean(t);
    w = 2*pi*f;
    A = [ones(size(t)), tau, cos(w*t), sin(w*t)];
    coeff = A\x;
    X = coeff(3) - 1j*coeff(4);   % complex phasor from cos/sin coefficients
end

function e = errMat(X, Xref)
    e = 100*norm(X-Xref,'fro') / max(norm(Xref,'fro'), eps);
end

function e = errComplex(x, xref)
    e = 100*abs(x-xref) / max(abs(xref), eps);
end

function printYmatrix(Y, names, pos)
    for e = 1:4
        r = pos(e,1); c = pos(e,2);
        fprintf('%s = %.9f < %.6f deg\n', names{e}, abs(Y(r,c)), rad2deg(angle(Y(r,c))));
    end
end