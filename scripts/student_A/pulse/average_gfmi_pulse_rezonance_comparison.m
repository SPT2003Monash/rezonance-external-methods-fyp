% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares accepted pulse-based average-GFMI admittance estimates against Rezonance results.

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
Ybase = 2e6/0.69e3^2;   % S (Sbase = 2 MVA, Vbase = 0.69 kV)

nrm = @(s) lower(regexprep(string(s),'[^a-zA-Z0-9]',''));
col = @(T,n) T.(T.Properties.VariableNames{find(ismember(nrm(T.Properties.VariableNames),nrm(n)),1)});
num = @(x) str2double(lower(string(x)));

pulseResultFile = fullfile(dataFolder,'AVG_pulse_admittance_accepted_24s.csv');
rezonanceFile = fullfile(dataFolder,'testing_rezo_0p5.csv');
if ~isfile(pulseResultFile)
    error('Pulse result file not found:\n%s', pulseResultFile);
end
if ~isfile(rezonanceFile)
    error('Rezonance file not found:\n%s', rezonanceFile);
end
P = readtable(pulseResultFile,'VariableNamingRule','preserve');
R = readtable(rezonanceFile,'VariableNamingRule','preserve');
% Frequency ranges
fp = num(col(P,{'f','Freq_Hz','Frequency','Frequency_Hz','freq','Hz'}));
fr = num(col(R,{'f','Frequency','Frequency_Hz','freq','Hz'}));
P = P(fp>=0.5 & fp<=5.3,:);   fp = fp(fp>=0.5 & fp<=5.3);
R = R(fr>=0.5 & fr<=5.0,:);   fr = fr(fr>=0.5 & fr<=5.0);
[fp,i] = sort(fp);  P = P(i,:);

% Insert NaNs at gaps (rejected bins near 4 Hz) so no false line is drawn
g = find(diff(fp) > 0.075);
[fx,ix] = sort([fp; fp(g)+0.01]);

pn = {'Ydd','Yqd','Ydq','Yqq'};
rn = {'dd','qd','dq','qq'};

figure('Color','w','Position',[80 80 1450 820]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for k = 1:4
    m = num(col(P,{['Mag_' pn{k}],[pn{k} '_mag'],['Mag' pn{k}]}));
    y = [m; nan(size(g))];
    z = abs(Ybase*num(col(R,{rn{k},pn{k}})));   % pu -> S

    nexttile;
    plot(fx,y(ix),'b.-','LineWidth',1.1,'MarkerSize',8); hold on;
    plot(fr,z,'rs','LineWidth',1.2,'MarkerSize',6,'MarkerFaceColor','w');
    grid on; xlim([0.5 5.3]);
    xlabel('Frequency (Hz)');
    ylabel(sprintf('|Y_{%s}| (S)',rn{k}));
    title(sprintf('Y_{%s}',rn{k}));
    if k==1, legend('Pulse qualified bins','Rezonance scan','Location','best'); end
end
sgtitle('Pulse versus Rezonance dq-Admittance Magnitude');
outputFigure = fullfile(figureFolder, ...
    'average_gfmi_pulse_vs_rezonance_magnitude.png');

exportgraphics(gcf, outputFigure, 'Resolution', 300);

fprintf('\nSaved comparison figure to:\n%s\n', outputFigure);