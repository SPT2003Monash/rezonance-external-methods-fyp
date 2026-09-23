% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Average-GFMI Richardson extrapolation analysis

clear; clc; close all;

%% Settings

fPert = 5;                   % 1, 5 or 10 Hz
operatingPoint = 'OP2';      % 'OP1' or 'OP2'

tStart = 6;
tEnd = 10;
thetaWindow = [8 10];
f0 = 50;

%% Paths

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

dataRoot = fullfile(repoRoot,'data','student_a','average_gfmi');
resultRoot = fullfile(dataRoot,'results');

if strcmpi(operatingPoint,'OP1')
    dataFolder = fullfile(dataRoot,'op1');
    prefix = 'AVG_GFMI';
    opLabel = 'OP1';
    pref = 1.0;
elseif strcmpi(operatingPoint,'OP2')
    dataFolder = fullfile(dataRoot,'op2');
    prefix = 'AVG_GFMI_OP2';
    opLabel = 'OP2';
    pref = 0.6;
else
    error('operatingPoint must be OP1 or OP2.');
end

resultFolder = fullfile(resultRoot,lower(opLabel));

if ~isfolder(resultFolder)
    mkdir(resultFolder);
end

freqText = sprintf('%gHz',fPert);

outputFile = fullfile(resultFolder,...
    sprintf('%s_%s_FINAL_results.csv',prefix,freqText));

fprintf('\n=== AVERAGE GFMI %s, Pref = %.1f pu, %g Hz ===\n',...
    opLabel,pref,fPert);

%% Load CSV files

axesTag = {'d','q'};
stepTag = {'h50','h25','h12p5','ref'};

fileStep.h50 = '50';
fileStep.h25 = '25';
fileStep.h12p5 = '12p5';
fileStep.ref = '0p1';

for a = 1:2
    ax = axesTag{a};

    for k = 1:numel(stepTag)
        step = stepTag{k};

        fileName = sprintf('%s_%s_%spert_%sus.csv',...
            prefix,freqText,ax,fileStep.(step));

        filePath = fullfile(dataFolder,fileName);

        fprintf('Loading: %s\n',filePath);
        D.(ax).(step) = loadPscadCsv(filePath);
    end
end

%% Check time grids

t = D.d.h50.Time;

for a = 1:2
    ax = axesTag{a};

    for k = 1:numel(stepTag)
        step = stepTag{k};
        tk = D.(ax).(step).Time;

        if numel(tk) ~= numel(t) || max(abs(tk-t)) > 1e-12
            error('%s %s file has a different time grid.',ax,step);
        end
    end
end

fprintf('Logging step = %.3f us\n',median(diff(t))*1e6);

%% Extract abc waveforms and apply Richardson

for a = 1:2
    ax = axesTag{a};

    for k = 1:numel(stepTag)
        step = stepTag{k};

        Vabc.(ax).(step) = getAbc(D.(ax).(step),'Vpoc');
        Iabc.(ax).(step) = getAbc(D.(ax).(step),'Ipoc');
    end

    Vre = richardsonSet(Vabc.(ax).h50,...
                        Vabc.(ax).h25,...
                        Vabc.(ax).h12p5);

    Ire = richardsonSet(Iabc.(ax).h50,...
                        Iabc.(ax).h25,...
                        Iabc.(ax).h12p5);

    Vabc.(ax).RE1 = Vre.RE1;
    Vabc.(ax).RE2 = Vre.RE2;
    Vabc.(ax).mixed = Vre.mixed;

    Iabc.(ax).RE1 = Ire.RE1;
    Iabc.(ax).RE2 = Ire.RE2;
    Iabc.(ax).mixed = Ire.mixed;
end

%% Operating-point check

idx = t >= tStart & t < tEnd;
tUse = t(idx);

fprintf('\n=== OPERATING-POINT CHECK ===\n');

for a = 1:2
    ax = axesTag{a};

    Vref = phaseRms(Vabc.(ax).ref(idx,:));
    Iref = phaseRms(Iabc.(ax).ref(idx,:));

    fprintf('\n%s perturbation\n',upper(ax));

    for k = 1:3
        step = stepTag{k};

        Vrms = phaseRms(Vabc.(ax).(step)(idx,:));
        Irms = phaseRms(Iabc.(ax).(step)(idx,:));

        dV = 100*abs(Vrms-Vref)/abs(Vref);
        dI = 100*abs(Irms-Iref)/abs(Iref);

        fprintf('%-6s dV = %.6f%%   dI = %.6f%%\n',...
            upper(step),dV,dI);
    end
end

%% Common dq angle

thetaIdx = t >= thetaWindow(1) & t < thetaWindow(2);

thetaD = estimateTheta(Vabc.d.h25,t,thetaIdx,2*pi*f0);
thetaQ = estimateTheta(Vabc.q.h25,t,thetaIdx,2*pi*f0);

thetaOffset = angle(exp(1j*thetaD)+exp(1j*thetaQ));
theta = 2*pi*f0*t + thetaOffset;

fprintf('\nCommon dq offset = %.6f deg\n',rad2deg(thetaOffset));

%% dq transformation and phasor extraction

cases = {'h50','h25','h12p5','RE1','RE2','mixed','ref'};

for a = 1:2
    ax = axesTag{a};

    for k = 1:numel(cases)
        c = cases{k};

        [Vd,Vq] = abcToDq(Vabc.(ax).(c),theta);
        [Id,Iq] = abcToDq(Iabc.(ax).(c),theta);

        ph.(ax).(c).Vd = getPhasor(tUse,Vd(idx),fPert);
        ph.(ax).(c).Vq = getPhasor(tUse,Vq(idx),fPert);
        ph.(ax).(c).Id = getPhasor(tUse,Id(idx),fPert);
        ph.(ax).(c).Iq = getPhasor(tUse,Iq(idx),fPert);
    end
end

%% Calculate dq-admittance matrices

fprintf('\n=== VOLTAGE MATRIX CONDITIONING ===\n');

for k = 1:numel(cases)
    c = cases{k};

    V.(c) = [ph.d.(c).Vd ph.q.(c).Vd;
             ph.d.(c).Vq ph.q.(c).Vq];

    I.(c) = [ph.d.(c).Id ph.q.(c).Id;
             ph.d.(c).Iq ph.q.(c).Iq];

    Y.(c) = I.(c)/V.(c);

    fprintf('%-7s cond(V) = %.6f\n',upper(c),cond(V.(c)));
end

%% Errors

names = ["Ydd";"Ydq";"Yqd";"Yqq"];
row = [1 1 2 2];
col = [1 2 1 2];

err50 = zeros(4,1);
err25 = zeros(4,1);
err12 = zeros(4,1);
errRE1 = zeros(4,1);
errRE2 = zeros(4,1);
errMixed = zeros(4,1);
p50_25 = zeros(4,1);
p25_12 = zeros(4,1);

for k = 1:4
    r = row(k);
    c = col(k);

    y50 = Y.h50(r,c);
    y25 = Y.h25(r,c);
    y12 = Y.h12p5(r,c);
    yref = Y.ref(r,c);

    err50(k) = errorPct(y50,yref);
    err25(k) = errorPct(y25,yref);
    err12(k) = errorPct(y12,yref);
    errRE1(k) = errorPct(Y.RE1(r,c),yref);
    errRE2(k) = errorPct(Y.RE2(r,c),yref);
    errMixed(k) = errorPct(Y.mixed(r,c),yref);

    p50_25(k) = log2(abs((y50-yref)/(y25-yref)));
    p25_12(k) = log2(abs((y25-yref)/(y12-yref)));
end

matrixErrors = [matrixError(Y.h50,Y.ref),...
                matrixError(Y.h25,Y.ref),...
                matrixError(Y.h12p5,Y.ref),...
                matrixError(Y.RE1,Y.ref),...
                matrixError(Y.RE2,Y.ref),...
                matrixError(Y.mixed,Y.ref)];

methodNames = {'50 us','25 us','12.5 us',...
               'p=1 Richardson','p=2 Richardson','Mixed Richardson'};

fprintf('\n=== FULL-MATRIX ERRORS ===\n');

for k = 1:numel(methodNames)
    fprintf('%-20s %.9f%%\n',methodNames{k},matrixErrors(k));
end

[bestError,bestIndex] = min(matrixErrors(4:6));
selectedCases = {'RE1','RE2','mixed'};
selectedNames = {'p=1 Richardson','p=2 Richardson',...
                 'Mixed h+h^2 Richardson'};

selectedCase = selectedCases{bestIndex};
selectedName = selectedNames{bestIndex};

fprintf('\nBest method: %s\n',selectedName);
fprintf('Full-matrix error: %.9f%%\n',bestError);

fprintf('\n=== DIRECT 0.1 us REFERENCE MATRIX ===\n');
printMatrix(Y.ref);

fprintf('\n=== SELECTED MATRIX: %s ===\n',selectedName);
printMatrix(Y.(selectedCase));

%% Save results

ResultsTable = table(names,err50,err25,err12,...
    errRE1,errRE2,errMixed,p50_25,p25_12,...
    'VariableNames',{'Element','Error_50us_pct',...
    'Error_25us_pct','Error_12p5us_pct',...
    'Error_RE1_pct','Error_RE2_pct',...
    'Error_Mixed_pct','Observed_p_50_to_25',...
    'Observed_p_25_to_12p5'});

writetable(ResultsTable,outputFile);

fprintf('\nSaved results to:\n%s\n',outputFile);

%% Plots

figure;
bar([err25 err12 errRE1 errRE2 errMixed]);
grid on;
set(gca,'XTickLabel',names);
ylabel('Error relative to 0.1 \mus reference (%)');
title(sprintf('%s, %g Hz: Element Errors',opLabel,fPert));
legend('25 \mus','12.5 \mus','p=1 RE','p=2 RE',...
       'Mixed RE','Location','best');

figure;
bar(matrixErrors(2:6));
grid on;
set(gca,'XTickLabel',methodNames(2:6));
ylabel('Full-matrix error (%)');
title(sprintf('%s, %g Hz: Full-Matrix Error',opLabel,fPert));

%% Local functions

function S = richardsonSet(h50,h25,h12)
S.RE1 = 2*h25-h50;
S.RE2 = (4*h25-h50)/3;
S.mixed = h50/3-2*h25+(8/3)*h12;
end

function D = loadPscadCsv(fileName)
if ~isfile(fileName)
    error('File not found:\n%s',fileName);
end

D = readtable(fileName,'VariableNamingRule','preserve');
D.Properties.VariableNames = strtrim(D.Properties.VariableNames);
end

function X = getAbc(D,prefix)
X = [D.([prefix '_A']),...
     D.([prefix '_B']),...
     D.([prefix '_C'])];
end

function value = phaseRms(x)
value = mean(sqrt(mean(x.^2,1)));
end

function [alpha,beta] = clarke(abc)
a = abc(:,1);
b = abc(:,2);
c = abc(:,3);

alpha = (2/3)*(a-0.5*b-0.5*c);
beta = (2/3)*(sqrt(3)/2)*(b-c);
end

function [d,q] = abcToDq(abc,theta)
[alpha,beta] = clarke(abc);

d = alpha.*cos(theta)+beta.*sin(theta);
q = -alpha.*sin(theta)+beta.*cos(theta);
end

function X = getPhasor(t,x,f)
tau = t-t(1);
w = 2*pi*f;

H = [cos(w*tau),sin(w*tau),ones(size(tau)),tau];
coef = H\x;

X = coef(1)-1j*coef(2);
end

function offset = estimateTheta(Vabc,t,idx,w0)
[alpha,beta] = clarke(Vabc);
spaceVector = alpha+1j*beta;

offset = angle(mean(spaceVector(idx).*exp(-1j*w0*t(idx))));
end

function e = errorPct(y,yref)
e = 100*abs(y-yref)/abs(yref);
end

function e = matrixError(Y,Yref)
e = 100*norm(Y-Yref,'fro')/norm(Yref,'fro');
end

function printMatrix(Y)
labels = {'Ydd','Ydq','Yqd','Yqq'};
positions = [1 1;1 2;2 1;2 2];

for k = 1:4
    r = positions(k,1);
    c = positions(k,2);

    fprintf('%s = %.9f < %.6f deg\n',...
        labels{k},abs(Y(r,c)),rad2deg(angle(Y(r,c))));
end
end