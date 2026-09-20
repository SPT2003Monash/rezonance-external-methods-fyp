% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares raw-waveform and envelope-based Richardson extrapolation for the switching GFMI model.

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
%% 1. SETTINGS --------------------------------------------------------------
freqList   = [1 5 10];

f0Hz      = 50;
pwmFreqHz = 5000;
lpfCutHz  = 100;

richOrder = 2;
stepRatio = 2;
richDenom = stepRatio^richOrder - 1;

thetaFitWindow  = [0.5 1.9];
phasorFitWindow = [3.0 10.0];

zoomFreqHz  = 10;      % which frequency's waveform to zoom in on
zoomAxis    = 'd';     % 'd' or 'q'
zoomSignal  = 'Id';    % 'Vd','Vq','Id','Iq'
zoomWindow  = [3.0 3.02];  % 20 ms window inside the perturbation

%% 2. RUN ALL THREE FREQUENCIES ----------------------------------------------

nF = numel(freqList);
rawErr = zeros(nF,3);   % columns: coarse, fine, RE
envErr = zeros(nF,3);

results = struct();

for i = 1:nF
    testFreqHz = freqList(i);

    [rawRes,envRes,t] = processFrequency(dataFolder,testFreqHz,f0Hz, ...
        pwmFreqHz,lpfCutHz,richDenom,thetaFitWindow,phasorFitWindow);

    rawErr(i,:) = [rawRes.errCoarse rawRes.errFine rawRes.errRE];
    envErr(i,:) = [envRes.errCoarse envRes.errFine envRes.errRE];

    results(i).freq   = testFreqHz;
    results(i).rawRes = rawRes;
    results(i).envRes = envRes;
    results(i).t      = t;

    fprintf('%2d Hz done. Raw RE = %.4f%%, Envelope RE = %.4f%%\n', ...
        testFreqHz,rawRes.errRE,envRes.errRE);
end

%% 3. TREND FIGURE — ERROR vs FREQUENCY, ALL LEVELS, RAW vs ENVELOPE --------

figure('Name','Error vs frequency: raw vs envelope');

semilogy(freqList,rawErr(:,1),'--o','Color',[0.85 0.33 0.10],'LineWidth',1.3); hold on;
semilogy(freqList,rawErr(:,2),'--s','Color',[0.85 0.33 0.10],'LineWidth',1.3);
semilogy(freqList,rawErr(:,3),'--^','Color',[0.85 0.33 0.10],'LineWidth',1.8);

semilogy(freqList,envErr(:,1),'-o','Color',[0.00 0.45 0.74],'LineWidth',1.3);
semilogy(freqList,envErr(:,2),'-s','Color',[0.00 0.45 0.74],'LineWidth',1.3);
semilogy(freqList,envErr(:,3),'-^','Color',[0.00 0.45 0.74],'LineWidth',1.8);

xlabel('Perturbation frequency (Hz)');
ylabel('Full-matrix admittance error (%, log scale)');
title('Switching GFMI: raw (dashed) vs envelope (solid) across all test cases');
legend({'Raw 12.5us','Raw 6.25us','Raw RE', ...
         'Envelope 12.5us','Envelope 6.25us','Envelope RE'}, ...
         'Location','northwest');
xticks(freqList);
grid on; box on;

%% 4. WAVEFORM ZOOM — WHAT RICHARDSON ACTUALLY SUBTRACTS --------------------

zi = find(freqList == zoomFreqHz,1);
if isempty(zi)
    error('zoomFreqHz must be one of freqList.');
end

t          = results(zi).t;
dqRaw      = results(zi).rawRes.dq;
dqEnv      = results(zi).envRes.dq;

zoomIdx = t >= zoomWindow(1) & t <= zoomWindow(2);

rawCoarse = dqRaw.(zoomAxis).coarse.(zoomSignal)(zoomIdx);
rawFine   = dqRaw.(zoomAxis).fine.(zoomSignal)(zoomIdx);
rawDiff   = rawFine - rawCoarse;   % this is the term RE divides by richDenom

envCoarse = dqEnv.(zoomAxis).coarse.(zoomSignal)(zoomIdx);
envFine   = dqEnv.(zoomAxis).fine.(zoomSignal)(zoomIdx);
envDiff   = envFine - envCoarse;

tz = t(zoomIdx)*1e3;   % ms for readability

figure('Name','Raw vs envelope: coarse-fine difference (what RE subtracts)');


plot(tz,rawDiff,'Color',[0.85 0.33 0.10],'LineWidth',1.1); hold on;
plot(tz,envDiff,'Color',[0.00 0.45 0.74],'LineWidth',1.3);
xlabel('Time (ms)');
ylabel(sprintf('%s fine - coarse',zoomSignal));
title(sprintf('%.0f Hz, %s axis: raw (orange) vs envelope (blue) difference signal', ...
    zoomFreqHz,zoomAxis));
legend('Raw','Envelope','Location','best');
grid on; box on;


rawP2P = max(rawDiff) - min(rawDiff);
envP2P = max(envDiff) - min(envDiff);

fprintf('\nRaw diff peak-to-peak = %.6g\n',rawP2P);
fprintf('Envelope diff peak-to-peak = %.6g\n',envP2P);
fprintf('Ratio (raw/envelope) = %.3f\n',rawP2P/envP2P);

%% ======================== LOCAL FUNCTIONS ===============================

function [rawRes,envRes,t] = processFrequency(dataFolder,testFreqHz,f0Hz, ...
    pwmFreqHz,lpfCutHz,richDenom,thetaFitWindow,phasorFitWindow)

    if ~ismember(testFreqHz,[1 5 10])
        error('testFreqHz must be 1, 5, or 10 Hz.');
    end

    fTag = sprintf('%dHz',testFreqHz);
    axesTag = {'d','q'};
    stepTag = {'coarse','fine','ref'};
    stepFileTag = struct('coarse','12p5us','fine','6p25us','ref','0p1us');

    for a = 1:2
        ax = axesTag{a};
        for s = 1:3
            st = stepTag{s};
            files.(ax).(st) = fullfile(dataFolder, ...
                sprintf('Switching_GFMI_%s_%spert_%s.csv',fTag,ax,stepFileTag.(st)));
        end
    end

    for a = 1:2
        ax = axesTag{a};
        for s = 1:3
            st = stepTag{s};
            data.(ax).(st) = loadPscadCsv(files.(ax).(st));
        end
    end

    t  = data.d.coarse.time;
    dt = median(diff(t));
    fs = 1/dt;

    for a = 1:2
        ax = axesTag{a};
        for s = 1:3
            st = stepTag{s};
            tk = data.(ax).(st).time;
            if length(tk) ~= length(t) || max(abs(tk-t)) > 1e-12
                error('%s.%s is not time aligned for %d Hz.',ax,st,testFreqHz);
            end
        end
    end

    Nma = round((1/pwmFreqHz)/dt);
    [lpfB,lpfA] = butter(4,lpfCutHz/(fs/2));

    for a = 1:2
        ax = axesTag{a};
        for s = 1:3
            st = stepTag{s};
            data.(ax).(st).Venv = extractEnvelope(data.(ax).(st).Vabc,Nma,lpfB,lpfA);
            data.(ax).(st).Ienv = extractEnvelope(data.(ax).(st).Iabc,Nma,lpfB,lpfA);
        end
    end

    w0 = 2*pi*f0Hz;
    rawRes = runFullPipeline(data,'Vabc','Iabc',t,thetaFitWindow, ...
        phasorFitWindow,testFreqHz,richDenom,w0);
    envRes = runFullPipeline(data,'Venv','Ienv',t,thetaFitWindow, ...
        phasorFitWindow,testFreqHz,richDenom,w0);

end

function res = runFullPipeline(data,Vfield,Ifield,t,thetaFitWindow, ...
    phasorFitWindow,testFreqHz,richDenom,w0)
% One dq-admittance pipeline: common angle -> Park -> Richardson ->
% phasor fit -> 2x2 matrix -> errors. Called once for raw fields
% ('Vabc'/'Iabc') and once for envelope fields ('Venv'/'Ienv').
% Returns res.dq with the full time-domain traces per level, so callers
% can inspect the actual waveforms (not just the final % error).

    axesTag  = {'d','q'};
    stepTag  = {'coarse','fine','ref'};
    levels   = {'coarse','fine','RE','ref'};
    sigs     = {'Vd','Vq','Id','Iq'};

    thetaIdx = t >= thetaFitWindow(1) & t <= thetaFitWindow(2);
    fitIdx   = t >= phasorFitWindow(1) & t <  phasorFitWindow(2);
    tFit     = t(fitIdx);

    thetaD = estimateThetaOffset(data.d.fine.(Vfield),t,thetaIdx,w0);
    thetaQ = estimateThetaOffset(data.q.fine.(Vfield),t,thetaIdx,w0);
    thetaOffset = angle(exp(1j*thetaD) + exp(1j*thetaQ));
    theta = w0*t + thetaOffset;

    dq = struct();
    for a = 1:2
        ax = axesTag{a};
        for s = 1:3
            st = stepTag{s};
            [Vd,Vq,Id,Iq] = parkTransform( ...
                data.(ax).(st).(Vfield),data.(ax).(st).(Ifield),theta);
            dq.(ax).(st) = struct('Vd',Vd,'Vq',Vq,'Id',Id,'Iq',Iq);
        end
        for k = 1:numel(sigs)
            sig = sigs{k};
            dq.(ax).RE.(sig) = dq.(ax).fine.(sig) + ...
                (dq.(ax).fine.(sig) - dq.(ax).coarse.(sig))/richDenom;
        end
    end

    Y = struct();
    for L = 1:numel(levels)
        lv = levels{L};
        ph = struct();
        for a = 1:2
            ax = axesTag{a};
            for k = 1:numel(sigs)
                sig = sigs{k};
                ph.(ax).(sig) = fitPhasorLS( ...
                    tFit,dq.(ax).(lv).(sig)(fitIdx),testFreqHz);
            end
        end
        V = [ph.d.Vd, ph.q.Vd; ph.d.Vq, ph.q.Vq];
        I = [ph.d.Id, ph.q.Id; ph.d.Iq, ph.q.Iq];
        Y.(lv) = I / V;
    end

    res.theta     = thetaOffset;
    res.Y         = Y;
    res.dq        = dq;
    res.errCoarse = matrixErrorPct(Y.coarse,Y.ref);
    res.errFine   = matrixErrorPct(Y.fine,Y.ref);
    res.errRE     = matrixErrorPct(Y.RE,Y.ref);

    elementName = {'Ydd','Ydq';'Yqd','Yqq'};
    res.elemErr  = zeros(4,1);
    res.elemName = cell(4,1);
    n = 1;
    for row = 1:2
        for col = 1:2
            res.elemErr(n)  = complexErrorPct(Y.RE(row,col),Y.ref(row,col));
            res.elemName{n} = elementName{row,col};
            n = n + 1;
        end
    end
end

function S = loadPscadCsv(file)

    if ~isfile(file)
        error('File not found:\n%s',file);
    end

    T = readtable(file);
    T.Properties.VariableNames = strtrim(T.Properties.VariableNames);

    required = {'Time','Vpoc_A','Vpoc_B','Vpoc_C','Ipoc_A','Ipoc_B','Ipoc_C'};
    for k = 1:numel(required)
        if ~ismember(required{k},T.Properties.VariableNames)
            error('Missing %s in %s',required{k},file);
        end
    end

    S.time = T.Time;
    S.Vabc = [T.Vpoc_A, T.Vpoc_B, T.Vpoc_C];
    S.Iabc = [T.Ipoc_A, T.Ipoc_B, T.Ipoc_C];

end

function Xenv = extractEnvelope(Xabc,N,lpfB,lpfA)

    Xenv = zeros(size(Xabc));
    for k = 1:3
        movAvg = movmean(Xabc(:,k),N,'Endpoints','shrink');
        Xenv(:,k) = filtfilt(lpfB,lpfA,movAvg);
    end

end

function offset = estimateThetaOffset(Vabc,t,idx,w0)

    [alpha,beta] = clarkeTransform(Vabc(:,1),Vabc(:,2),Vabc(:,3));
    spaceVec = alpha + 1j*beta;
    offset = angle(mean(spaceVec(idx).*exp(-1j*w0*t(idx))));

end

function [Vd,Vq,Id,Iq] = parkTransform(Vabc,Iabc,theta)

    [Va,Vb] = clarkeTransform(Vabc(:,1),Vabc(:,2),Vabc(:,3));
    [Ia,Ib] = clarkeTransform(Iabc(:,1),Iabc(:,2),Iabc(:,3));

    Vd =  Va.*cos(theta) + Vb.*sin(theta);
    Vq = -Va.*sin(theta) + Vb.*cos(theta);
    Id =  Ia.*cos(theta) + Ib.*sin(theta);
    Iq = -Ia.*sin(theta) + Ib.*cos(theta);

end

function [alpha,beta] = clarkeTransform(a,b,c)

    alpha = sqrt(2/3)*(a - 0.5*b - 0.5*c);
    beta  = sqrt(2/3)*((sqrt(3)/2)*b - (sqrt(3)/2)*c);

end

function X = fitPhasorLS(t,x,f)

    t = t(:); x = x(:);
    tau = t - mean(t);
    w = 2*pi*f;

    A = [ones(size(t)), tau, cos(w*t), sin(w*t)];
    coeff = A\x;
    X = coeff(3) - 1j*coeff(4);

end

function err = matrixErrorPct(X,Xref)

    err = 100*norm(X-Xref,'fro') / max(norm(Xref,'fro'),eps);

end

function err = complexErrorPct(x,xref)

    err = 100*abs(x-xref) / max(abs(xref),eps);

end
