% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Evaluates first-order, second-order and mixed Richardson methods for the average GFMI model.

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
%% SETTINGS ------------------------------------------------------------

fPert = 10;                         % 1, 5 or 10 Hz
operatingPoint = '';            % Use '' for OP1, or 'OP2' for OP2

tStart = 3.0;
tEnd = 10.0;
thetaFitWindow = [0.5 1.9];
f0Hz = 50;
currentSign = +1;

if ~ismember(fPert,[1 5 10])
    error('fPert must be 1, 5, or 10 Hz.');
end

if isempty(operatingPoint)
    filePrefix = 'AVG_GFMI';
    opLabel = 'OP1';
    prefPu = 1.0;
    dataFolder = op1Folder;
elseif strcmpi(operatingPoint,'OP2')
    filePrefix = 'AVG_GFMI_OP2';
    opLabel = 'OP2';
    prefPu = 0.6;
    dataFolder = op2Folder;
else
    error('operatingPoint must be empty '''' for OP1 or ''OP2''.');
end
resultFolder = fullfile(resultRoot, lower(opLabel));

if ~isfolder(resultFolder)
    mkdir(resultFolder);
end

freqText = sprintf('%gHz',fPert);

outputFile = fullfile(resultFolder, sprintf('%s_%s_FINAL_RE1_RE2_Mixed_ref2_results.csv',filePrefix, freqText));

fprintf('\n=== VALIDATED AVERAGE GFMI %s, Pref = %.1f pu, %g Hz ===\n',opLabel,prefPu,fPert);
fprintf('Richardson inputs: 50/25/12.5 us | Held-out reference: 2 us\n');
fprintf('Theta fit: %.1f-%.1f s | Analysis window: %.1f-%.1f s\n',thetaFitWindow(1),thetaFitWindow(2),tStart,tEnd);

%% LOAD DATA + TIME-ALIGNMENT CHECK -------------------------------------

axesTag = {'d','q'};
stepTag = {'h50','h25','h12p5','ref'};
stepFileTag = struct('h50','50','h25','25','h12p5','12p5','ref','2');
fprintf('\n--- Loading files ---\n');

for a = 1:2
    for k = 1:numel(stepTag)

        ax = axesTag{a};
        s = stepTag{k};

        file = fullfile(dataFolder,sprintf('%s_%s_%spert_%sus.csv',filePrefix,freqText,ax,stepFileTag.(s)));
        fprintf('%s\n',file);
        D.(ax).(s) = loadPscadCsv(file);
    end
end

t = D.d.h50.Time;
for a = 1:2
    for k = 1:numel(stepTag)
        ax = axesTag{a}; s = stepTag{k};
        tk = D.(ax).(s).Time;
        if numel(tk) ~= numel(t) || max(abs(tk-t)) > 1e-12
            error('%s.%s is not time-aligned with the reference grid.',ax,s);
        end
    end
end
dtLog = median(diff(t));
fprintf('\nAll 8 files share the same time grid. Logging step = %.6f us (%.3f kHz)\n', dtLog*1e6, 1/dtLog/1000);

processingTimer = tic;

%% EXTRACT abc + RICHARDSON EXTRAPOLATION --------------------------------

for a = 1:2
    ax = axesTag{a};
    for k = 1:numel(stepTag)
        s = stepTag{k};
        Vabc.(ax).(s) = getAbc(D.(ax).(s),'Vpoc');
        Iabc.(ax).(s) = currentSign * getAbc(D.(ax).(s),'Ipoc');
    end

    reV = richardsonSet(Vabc.(ax).h50, Vabc.(ax).h25, Vabc.(ax).h12p5);
    reI = richardsonSet(Iabc.(ax).h50, Iabc.(ax).h25, Iabc.(ax).h12p5);
    Vabc.(ax).RE1 = reV.RE1; Vabc.(ax).RE2 = reV.RE2; Vabc.(ax).mixed = reV.mixed;
    Iabc.(ax).RE1 = reI.RE1; Iabc.(ax).RE2 = reI.RE2; Iabc.(ax).mixed = reI.mixed;
end

idx = t >= tStart & t < tEnd;
tUse = t(idx);
if nnz(idx) < 100; error('Too few samples in analysis window.'); end

%% OPERATING-POINT RMS CHECK ---------------------------------------------

fprintf('\n=== OPERATING-POINT RMS CHECK ===\n');
for a = 1:2
    ax = axesTag{a};
    fprintf('\n%s-axis perturbation\n',upper(ax));
    for k = 1:numel(stepTag)
        s = stepTag{k};
        Vrms = mean(sqrt(mean(Vabc.(ax).(s)(idx,:).^2,1)));
        Irms = mean(sqrt(mean(Iabc.(ax).(s)(idx,:).^2,1)));
        fprintf('%-7s : V = %.9f   I = %.9f\n', upper(s), Vrms, Irms);
    end
end

fprintf('\n--- Operating-point consistency vs 2 us ---\n');
for a = 1:2
    ax = axesTag{a};
    VrefRMS = mean(sqrt(mean(Vabc.(ax).ref(idx,:).^2,1)));
    IrefRMS = mean(sqrt(mean(Iabc.(ax).ref(idx,:).^2,1)));

    for k = 1:3   % h50, h25, h12p5 -- exclude ref itself
        s = stepTag{k};
        Vrms = mean(sqrt(mean(Vabc.(ax).(s)(idx,:).^2,1)));
        Irms = mean(sqrt(mean(Iabc.(ax).(s)(idx,:).^2,1)));
        dV = 100*abs(Vrms-VrefRMS)/max(abs(VrefRMS),eps);
        dI = 100*abs(Irms-IrefRMS)/max(abs(IrefRMS),eps);
        fprintf('%s %-6s : dV = %.6f %%   dI = %.6f %%\n', upper(ax), upper(s), dV, dI);
        if dV > 1 || dI > 1
            warning('%s %s operating point differs by more than 1%%.', upper(ax), upper(s));
        end
    end
end

%% COMMON dq ANGLE (from pre-perturbation 25 us runs only) ---------------

thetaIdx = t >= thetaFitWindow(1) & t <= thetaFitWindow(2);
thetaOffsetD = estimateThetaOffset(Vabc.d.h25, t, thetaIdx, 2*pi*f0Hz);
thetaOffsetQ = estimateThetaOffset(Vabc.q.h25, t, thetaIdx, 2*pi*f0Hz);
thetaOffset = angle(exp(1j*thetaOffsetD) + exp(1j*thetaOffsetQ));   % circular mean
theta = 2*pi*f0Hz*t + thetaOffset;

fprintf('\n=== COMMON dq ANGLE ===\n');
fprintf('25us d offset = %.9f rad (%.6f deg)\n', thetaOffsetD, rad2deg(thetaOffsetD));
fprintf('25us q offset = %.9f rad (%.6f deg)\n', thetaOffsetQ, rad2deg(thetaOffsetQ));
fprintf('Common offset = %.9f rad (%.6f deg)\n', thetaOffset, rad2deg(thetaOffset));

% Diagnostic only: confirms the 25us-based angle is close to what 2us gives
thetaRefD = estimateThetaOffset(Vabc.d.ref, t, thetaIdx, 2*pi*f0Hz);
thetaRefQ = estimateThetaOffset(Vabc.q.ref, t, thetaIdx, 2*pi*f0Hz);
fprintf('\n--- Angle diagnostic (25 us vs 2 us) ---\n');
fprintf('d difference = %.6f deg\n', wrappedAngleDifference(thetaRefD,thetaOffsetD));
fprintf('q difference = %.6f deg\n', wrappedAngleDifference(thetaRefQ,thetaOffsetQ));

%% dq TRANSFORM + PHASOR EXTRACTION ---------------------------------------

allCases = {'h50','h25','h12p5','RE1','RE2','mixed','ref'};
for a = 1:2
    ax = axesTag{a};
    for c = 1:numel(allCases)
        s = allCases{c};
        [Vd,Vq] = abcToDq(Vabc.(ax).(s), theta);
        [Id,Iq] = abcToDq(Iabc.(ax).(s), theta);
        ph.(ax).(s).Vd = estimatePhasor(tUse, Vd(idx), fPert);
        ph.(ax).(s).Vq = estimatePhasor(tUse, Vq(idx), fPert);
        ph.(ax).(s).Id = estimatePhasor(tUse, Id(idx), fPert);
        ph.(ax).(s).Iq = estimatePhasor(tUse, Iq(idx), fPert);
    end
end

%% FORM V, I, Y MATRICES ----------------------------------------------------

fprintf('\n=== VOLTAGE MATRIX CONDITIONING ===\n');
for c = 1:numel(allCases)
    s = allCases{c};
    Vmat.(s) = [ph.d.(s).Vd, ph.q.(s).Vd; ph.d.(s).Vq, ph.q.(s).Vq];
    Imat.(s) = [ph.d.(s).Id, ph.q.(s).Id; ph.d.(s).Iq, ph.q.(s).Iq];
    Y.(s) = Imat.(s) / Vmat.(s);   % right division -- avoids inv(V)
    fprintf('%-8s cond(V) = %.6f\n', upper(s), cond(Vmat.(s)));
end

%% ELEMENTWISE ERROR + CONVERGENCE -------------------------------------------

names = ["Ydd";"Ydq";"Yqd";"Yqq"];
rows = [1;1;2;2]; cols = [1;2;1;2];
nElem = 4;
err50=zeros(nElem,1); err25=zeros(nElem,1); err12=zeros(nElem,1);
errRE1=zeros(nElem,1); errRE2=zeros(nElem,1); errMix=zeros(nElem,1);
p50_25=zeros(nElem,1); p25_12=zeros(nElem,1);
angle50_25=zeros(nElem,1); angle25_12=zeros(nElem,1);

fprintf('\n=== ELEMENT RESULTS AT %g Hz (reference = direct 2 us) ===\n',fPert);
for k = 1:nElem
    rr = rows(k); cc = cols(k);
    y50=Y.h50(rr,cc); y25=Y.h25(rr,cc); y12=Y.h12p5(rr,cc);
    yRE1=Y.RE1(rr,cc); yRE2=Y.RE2(rr,cc); yMix=Y.mixed(rr,cc); yRef=Y.ref(rr,cc);

    err50(k)=errorPct(y50,yRef); err25(k)=errorPct(y25,yRef); err12(k)=errorPct(y12,yRef);
    errRE1(k)=errorPct(yRE1,yRef); errRE2(k)=errorPct(yRE2,yRef); errMix(k)=errorPct(yMix,yRef);

    ratio50_25 = (y50-yRef)/(y25-yRef);
    ratio25_12 = (y25-yRef)/(y12-yRef);
    p50_25(k) = log(abs(ratio50_25))/log(2);
    p25_12(k) = log(abs(ratio25_12))/log(2);
    angle50_25(k) = rad2deg(angle(ratio50_25));
    angle25_12(k) = rad2deg(angle(ratio25_12));

    fprintf('\n--- %s ---\n', names(k));
    fprintf('Direct: 50us=%.9f%%, 25us=%.9f%%, 12.5us=%.9f%%\n', err50(k), err25(k), err12(k));
    fprintf('Richardson: p=1=%.9f%%, p=2=%.9f%%, mixed=%.9f%%\n', errRE1(k), errRE2(k), errMix(k));
    fprintf('Observed p: 50->25=%.4f, 25->12.5=%.4f | angle: 50->25=%.3fdeg, 25->12.5=%.3fdeg\n', ...
        p50_25(k), p25_12(k), angle50_25(k), angle25_12(k));

    [bestElemErr,bestElemIdx] = min([errRE1(k), errRE2(k), errMix(k)]);
    elementREnames = {'p=1 Richardson','p=2 Richardson','Mixed h+h^2 Richardson'};
    fprintf('BEST RICHARDSON = %s (%.9f%%)\n', elementREnames{bestElemIdx}, bestElemErr);
end

%% FULL-MATRIX ERRORS AND COMPARISON -----------------------------------------

matrixErr50  = matrixError(Y.h50,  Y.ref);
matrixErr25  = matrixError(Y.h25,  Y.ref);
matrixErr12  = matrixError(Y.h12p5,Y.ref);
matrixErrRE1 = matrixError(Y.RE1,  Y.ref);
matrixErrRE2 = matrixError(Y.RE2,  Y.ref);
matrixErrMix = matrixError(Y.mixed,Y.ref);

fprintf('\n=== FULL-MATRIX COMPARISON AT %g Hz (reference = direct 2 us) ===\n',fPert);
fprintf('Direct PSCAD:        50us=%.9f%%, 25us=%.9f%%, 12.5us=%.9f%%\n', matrixErr50, matrixErr25, matrixErr12);
fprintf('Waveform Richardson: p=1=%.9f%%, p=2=%.9f%%, mixed=%.9f%%\n', matrixErrRE1, matrixErrRE2, matrixErrMix);

fprintf('\n=== IMPROVEMENT FACTORS ===\n');
fprintf('p=1 RE vs direct 25us     = %.3fx\n', matrixErr25/max(matrixErrRE1,eps));
fprintf('p=2 RE vs direct 25us     = %.3fx\n', matrixErr25/max(matrixErrRE2,eps));
fprintf('Mixed RE vs direct 12.5us = %.3fx\n', matrixErr12/max(matrixErrMix,eps));

%% SELECT BEST RICHARDSON METHOD ----------------------------------------------

REerrors = [matrixErrRE1, matrixErrRE2, matrixErrMix];
REnames = {'p=1 Richardson','p=2 Richardson','Mixed h+h^2 Richardson'};
REmatrices = {Y.RE1, Y.RE2, Y.mixed};
[bestREerror,bestREidx] = min(REerrors);
selectedName = REnames{bestREidx};
Yselected = REmatrices{bestREidx};

fprintf('\n=== BEST RICHARDSON METHOD ===\n');
fprintf('Best method = %s, error = %.9f%%\n', selectedName, bestREerror);

fprintf('\n=== MIXED-RICHARDSON CHECK ===\n');
fprintf('p=1=%.9f%%, p=2=%.9f%%, mixed=%.9f%%\n', matrixErrRE1, matrixErrRE2, matrixErrMix);
fprintf('Mixed beats p=1: %s | Mixed beats p=2: %s\n', ...
    yesNo(matrixErrMix < matrixErrRE1), yesNo(matrixErrMix < matrixErrRE2));

fprintf('\n=== DIRECT 2 us REFERENCE MATRIX ===\n'); printYmatrix(Y.ref);
fprintf('\n=== SELECTED RICHARDSON MATRIX (%s) ===\n', selectedName); printYmatrix(Yselected);

%% RESULTS TABLE ---------------------------------------------------------------

ResultsTable = table(names, err50,err25,err12, errRE1,errRE2,errMix, ...
    p50_25,p25_12, angle50_25,angle25_12, ...
    'VariableNames',{'Element','Error_50us_pct','Error_25us_pct','Error_12p5us_pct', ...
    'Error_RE1_pct','Error_RE2_pct','Error_Mixed_pct', ...
    'Observed_p_50_to_25','Observed_p_25_to_12p5', ...
    'ErrorAngle_50_to_25_deg','ErrorAngle_25_to_12p5_deg'});

disp(ResultsTable);
writetable(ResultsTable, outputFile);
fprintf('\nSaved results to:\n%s\n', outputFile);

matlabProcessingTime = toc(processingTimer);
fprintf('\nMATLAB post-processing time = %.6f s\n', matlabProcessingTime);

%% PLOTS -------------------------------------------------------------------

figure;
bar([err25, err12, errRE1, errRE2, errMix]);
grid on; set(gca,'XTickLabel',names);
ylabel('Complex error vs 2 us (%)');
title(sprintf('%g Hz: Richardson method comparison', fPert));
legend('25 us direct','12.5 us direct','p=1 RE','p=2 RE','Mixed h+h^2 RE','Location','best');

figure;
bar([matrixErr25, matrixErr12, matrixErrRE1, matrixErrRE2, matrixErrMix]);
grid on; set(gca,'XTickLabel',{'25 us','12.5 us','p=1 RE','p=2 RE','Mixed RE'});
ylabel('Full-matrix error vs 2 us (%)');
title(sprintf('%g Hz: Full-matrix Richardson comparison', fPert));

%% ============================ LOCAL FUNCTIONS ================================

function S = richardsonSet(h50,h25,h12p5)
    S.RE1   = 2*h25 - h50;                       % 1st-order, cancels O(h)
    S.RE2   = (4*h25 - h50)/3;                   % 2nd-order, cancels O(h^2)
    S.mixed = (1/3)*h50 - 2*h25 + (8/3)*h12p5;   % cancels O(h) and O(h^2)
end

function D = loadPscadCsv(fileName)
    if ~isfile(fileName); error('File not found:\n%s', fileName); end
    D = readtable(fileName, 'VariableNamingRule','preserve');
    D.Properties.VariableNames = strtrim(D.Properties.VariableNames);

    required = {'Time','Vpoc_A','Vpoc_B','Vpoc_C','Ipoc_A','Ipoc_B','Ipoc_C'};
    for k = 1:numel(required)
        if ~ismember(required{k}, D.Properties.VariableNames)
            error('Missing column "%s" in %s.', required{k}, fileName);
        end
    end
end

function X = getAbc(D,prefix)
    X = [D.([prefix '_A']), D.([prefix '_B']), D.([prefix '_C'])];
end

function [alpha,beta] = clarkeTransform(abc)
    a = abc(:,1); b = abc(:,2); c = abc(:,3);
    alpha = (2/3)*(a - 0.5*b - 0.5*c);            % amplitude-invariant Clarke transform
    beta  = (2/3)*((sqrt(3)/2)*b - (sqrt(3)/2)*c);
end

function [d,q] = abcToDq(abc,theta)
    [alpha,beta] = clarkeTransform(abc);
    d =  alpha.*cos(theta) + beta.*sin(theta);
    q = -alpha.*sin(theta) + beta.*cos(theta);
end

function X = estimatePhasor(t,x,f)
    tau = t - t(1);
    w = 2*pi*f;
    H = [cos(w*tau), sin(w*tau), ones(size(tau)), tau];   % + DC & trend terms
    beta = H\x;
    X = beta(1) - 1j*beta(2);   % complex phasor: x(t) = Re{X exp(jwt)}
end

function e = errorPct(y,yref)
    e = 100*abs(y-yref) / max(abs(yref), eps);
end

function e = matrixError(Y,Yref)
    e = 100*norm(Y-Yref,'fro') / max(norm(Yref,'fro'), eps);
end

function offset = estimateThetaOffset(Vabc,t,idx,w0)
    [alpha,beta] = clarkeTransform(Vabc);
    spaceVec = alpha + 1j*beta;
    offset = angle(mean(spaceVec(idx) .* exp(-1j*w0*t(idx))));
end

function ddeg = wrappedAngleDifference(a,b)
    ddeg = rad2deg(angle(exp(1j*(a-b))));
end

function txt = yesNo(condition)
    if condition; txt = 'YES'; else; txt = 'NO'; end
end

function printYmatrix(Y)
    elemNames = {'Ydd','Ydq','Yqd','Yqq'};
    positions = [1 1; 1 2; 2 1; 2 2];
    for k = 1:4
        r = positions(k,1); c = positions(k,2);
        fprintf('%s = %.9f < %.6f deg\n', elemNames{k}, abs(Y(r,c)), rad2deg(angle(Y(r,c))));
    end
end