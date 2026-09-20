% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares selected average-GFMI Richardson admittance matrices against Rezonance results.

clear; clc; close all;
%% Repository paths
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

averageDataRoot = fullfile(repoRoot, ...
    'data', 'student_a', 'average_gfmi');

op1Folder = fullfile(averageDataRoot, 'op1');

resultFolder = fullfile(averageDataRoot, ...
    'results', 'rezonance_comparison');

figureFolder = fullfile(repoRoot, ...
    'figures', 'student_a', 'average_gfmi');

if ~isfolder(op1Folder)
    error('OP1 data folder not found:\n%s', op1Folder);
end

if ~isfolder(resultFolder)
    mkdir(resultFolder);
end

if ~isfolder(figureFolder)
    mkdir(figureFolder);
end

%% BASE ADMITTANCE ----------------------------------------------------------

Sbase = 2e6;    % VA
Vbase = 690;    % V line-line
Ybase = Sbase / Vbase^2;

fprintf('\n=== BASE ADMITTANCE ===\n');
fprintf('Sbase = %.3f MVA, Vbase = %.3f kV, Ybase = %.9f S\n', Sbase/1e6, Vbase/1000, Ybase);

%% REZONANCE DATA -------------------------------------------------------------

freq  = [1 5 10];
names = {'Ydd','Ydq','Yqd','Yqq'};
method = {'p=1 Richardson (50 + 25 us)', 'Mixed h+h^2 RE (50 + 25 + 12.5 us)', 'Mixed h+h^2 RE (50 + 25 + 12.5 us)'};

rezonanceFile = fullfile(op1Folder, 'rezonance_inverter.csv');
if ~isfile(rezonanceFile)
    error('Rezonance input file not found: %s', rezonanceFile);
end

RZraw = readtable(rezonanceFile, 'VariableNamingRule','preserve');
RZall = [parseComplexCol(RZraw.dd), parseComplexCol(RZraw.dq), parseComplexCol(RZraw.qd), parseComplexCol(RZraw.qq)];

[~,rowIdx] = ismember(freq, RZraw.f);
if any(rowIdx==0)
    error('Rezonance CSV is missing one of the requested frequencies: %s Hz', mat2str(freq(rowIdx==0)));
end
RZ = RZall(rowIdx,:);              % rows = 1/5/10 Hz, columns = Ydd/Ydq/Yqd/Yqq
RZ_mag = abs(RZ);
RZ_phase = rad2deg(angle(RZ));

%% SELECTED FINAL RICHARDSON RESULTS (Siemens, final corrected values) --------

RE_mag_S = [ 4.370935106,  6.520009587, 10.456836737, 8.290797808;
             6.348082724, 13.068358501,  9.392693991, 9.657905288;
             8.000908668, 13.788420374, 10.557171114, 9.011811423];

RE_phase_original = [-172.570424, -122.807407, -2.390393,  -54.329041;
                      -159.062767, -175.199092, -7.712001, -151.414814;
                      -154.878902,  169.596835, -9.409932, -159.237520];

RE_S = RE_mag_S .* exp(1j*deg2rad(RE_phase_original));

%% CONVERT TO pu + CURRENT-DIRECTION CORRECTION ------------------------------

RE_corrected = -(RE_S / Ybase);
RE_mag = abs(RE_corrected);
RE_phase = rad2deg(angle(RE_corrected));

%% DIFFERENCES ----------------------------------------------------------------

magDifferencePct = 100 * abs(RE_mag - RZ_mag) ./ RZ_mag;
phaseDifferenceDeg = rad2deg(angle(exp(1j*deg2rad(RE_phase - RZ_phase))));
complexErrorPct = 100 * abs(RE_corrected - RZ) ./ abs(RZ);

%% PRINT ELEMENT RESULTS -------------------------------------------------------

fprintf('\n=== REZONANCE vs SELECTED RICHARDSON ===\n');
for f = 1:numel(freq)
    fprintf('\n--- %g Hz (%s) ---\n', freq(f), method{f});
    for k = 1:4
        fprintf('%s: RZ = %.5f<%.3fdeg | RE = %.5f<%.3fdeg | dMag=%.3f%% dPhase=%.3fdeg complexErr=%.3f%%\n', ...
            names{k}, RZ_mag(f,k), RZ_phase(f,k), RE_mag(f,k), RE_phase(f,k), ...
            magDifferencePct(f,k), phaseDifferenceDeg(f,k), complexErrorPct(f,k));
    end
end

%% FULL 2x2 MATRIX ERROR (Frobenius norm) ---------------------------------------

matrixErrorPct = zeros(3,1);
fprintf('\n=== FULL-MATRIX ERROR vs REZONANCE ===\n');
for f = 1:3
    Yrz = reshape(RZ(f,:),2,2).';
    Yre = reshape(RE_corrected(f,:),2,2).';
    matrixErrorPct(f) = 100*norm(Yre-Yrz,'fro') / norm(Yrz,'fro');
    fprintf('%g Hz = %.6f %%\n', freq(f), matrixErrorPct(f));
end

%% SAVE RESULTS TABLES ----------------------------------------------------------

Frequency = repelem(freq',4);
Element = repmat(string(names)',3,1);

ComparisonTable = table(Frequency, Element, ...
    reshape(RZ_mag.',[],1), reshape(RZ_phase.',[],1), ...
    reshape(RE_mag.',[],1), reshape(RE_phase.',[],1), ...
    reshape(magDifferencePct.',[],1), reshape(phaseDifferenceDeg.',[],1), reshape(complexErrorPct.',[],1), ...
    'VariableNames',{'Frequency','Element','Rezonance_Magnitude_pu','Rezonance_Phase_deg', ...
    'Richardson_Magnitude_pu','Richardson_Phase_deg','Magnitude_Difference_pct', ...
    'Phase_Difference_deg','Complex_Error_pct'});

fprintf('\n=== FINAL COMPARISON TABLE ===\n\n');
disp(ComparisonTable);
comparisonFile = fullfile(resultFolder, 'Average_GFMI_Rezonance_vs_Richardson.csv');
writetable(ComparisonTable, comparisonFile);

MatrixErrorTable = table(freq', matrixErrorPct, 'VariableNames',{'Frequency_Hz','FullMatrix_Error_pct'});
matrixErrorFile = fullfile(resultFolder,'Average_GFMI_Rezonance_vs_Richardson_MatrixError.csv');
writetable(MatrixErrorTable, matrixErrorFile);

fprintf('\nSaved:\n%s\n%s\n', comparisonFile, matrixErrorFile);

%% PLOTS: MAGNITUDE AND PHASE ----------------------------------------------------

plotComparison(RZ_mag, RE_mag, names, freq,'Admittance magnitude (pu)','Rezonance vs Richardson - Magnitude',fullfile(figureFolder,'average_gfmi_rezonance_vs_richardson_magnitude.png'));
plotComparison(RZ_phase, RE_phase, names, freq,'Phase (deg)','Rezonance vs Richardson - Phase',fullfile(figureFolder,'average_gfmi_rezonance_vs_richardson_phase.png'));
%% ============================ LOCAL FUNCTIONS ================================
function plotComparison(rz, re, names, freq, yLabelText, figTitle, outputFile)
    figure('Name', figTitle);
    for f = 1:numel(freq)
        subplot(1,numel(freq),f);
        bar([rz(f,:)', re(f,:)']);
        grid on; set(gca,'XTickLabel',names);
        ylabel(yLabelText);
        title(sprintf('%g Hz', freq(f)));
        if f == 1; legend('Rezonance','Selected Richardson','Location','best'); end
    end
    sgtitle(figTitle);
    exportgraphics(gcf, outputFile, 'Resolution', 300);
end

function z = parseComplexCol(col)
    % Parses a table column of Python-style complex strings, e.g. '0.87+0.26j',
    % into a MATLAB complex double column vector.
    if isnumeric(col); z = double(col); return; end
    col = string(col);
    z = zeros(numel(col),1);
    for k = 1:numel(col)
        z(k) = parseComplexStr(col(k));
    end
end

function z = parseComplexStr(s)
    s = strtrim(char(s));
    tok = regexp(s, '^([+-]?[\d.]+(?:[eE][+-]?\d+)?)([+-][\d.]+(?:[eE][+-]?\d+)?)j$', 'tokens', 'once');
    if isempty(tok); error('Could not parse complex value: %s', s); end
    z = str2double(tok{1}) + 1i*str2double(tok{2});
end