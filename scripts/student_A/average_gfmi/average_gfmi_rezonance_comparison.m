% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares average-GFMI Richardson and Rezonance results.

clear; clc; close all;

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(scriptDir)));

dataFile = fullfile(repoRoot,'data','student_a','average_gfmi','op1','rezonance_inverter.csv');
resultDir = fullfile(repoRoot,'data','student_a','average_gfmi','results','rezonance_comparison');
figureDir = fullfile(repoRoot,'figures','student_a','average_gfmi');

if ~isfolder(resultDir); mkdir(resultDir); end
if ~isfolder(figureDir); mkdir(figureDir); end

%% Rezonance results
freq = [1 5 10];
names = {'Ydd','Ydq','Yqd','Yqq'};
Ybase = 2e6/690^2;

T = readtable(dataFile,'VariableNamingRule','preserve');
RZpu = [toComplex(T.dd),toComplex(T.dq),toComplex(T.qd),toComplex(T.qq)];
[~,idx] = ismember(freq,T.f);
RZ = -Ybase*RZpu(idx,:);       % pu to Siemens and current-direction correction

%% Direct 0.1 us reference
refMag = [
 3.805087893  6.533549111 10.970783680 7.762708222
 6.336337503 12.241892992  9.094264464 9.455390447
 7.906670048 13.136877485 10.227463818 9.079387525];

refPhase = [
-164.836334 -118.939498 -4.758002  -61.323822
-157.951183 -174.244815 -8.813936 -153.790456
-154.421456  171.358587 -9.341142 -161.496981];

REF = refMag.*exp(1j*deg2rad(refPhase));

%% Selected mixed Richardson results
reMag = [
 3.805187663  6.533628624 10.970854219 7.762753468
 6.336246268 12.241927847  9.094203690 9.455275573
 7.906333252 13.136621340 10.227415340 9.079086044];

rePhase = [
-164.835375 -118.939577 -4.758411  -61.323817
-157.953834 -174.246412 -8.814431 -153.792448
-154.423517  171.356349 -9.342690 -161.498808];

RE = reMag.*exp(1j*deg2rad(rePhase));

%% Errors against direct 0.1 us reference
elementRE = 100*abs(RE-REF)./abs(REF);
elementRZ = 100*abs(RZ-REF)./abs(REF);

matrixRE = vecnorm(RE-REF,2,2)./vecnorm(REF,2,2)*100;
matrixRZ = vecnorm(RZ-REF,2,2)./vecnorm(REF,2,2)*100;

fprintf('\nFrequency   Mixed RE error   Rezonance error\n');
for k = 1:3
    fprintf('%5g Hz      %10.6f %%      %10.6f %%\n', ...
        freq(k),matrixRE(k),matrixRZ(k));
end

%% Save results
Results = table(freq',matrixRE,matrixRZ,...
    'VariableNames',{'Frequency_Hz','Mixed_Error_pct','Rezonance_Error_pct'});

writetable(Results,fullfile(resultDir,...
    'Average_GFMI_Rezonance_comparison.csv'));

%% Plots
figure;
for k = 1:3
    subplot(1,3,k);
    bar([abs(REF(k,:));abs(RE(k,:));abs(RZ(k,:))]');
    grid on; title(sprintf('%g Hz',freq(k)));
    set(gca,'XTickLabel',names);
    if k == 1; ylabel('Admittance (S)'); end
end
legend('Direct 0.1 \mus','Mixed Richardson','Rezonance');
sgtitle('Average GFMI admittance magnitude');
exportgraphics(gcf,fullfile(figureDir,...
    'average_gfmi_rezonance_magnitude.png'),'Resolution',300);

%% Phase comparison
refPhase = rad2deg(angle(REF));
rePhase  = rad2deg(angle(RE));
rzPhase  = rad2deg(angle(RZ));

figure;
for k = 1:3
    subplot(1,3,k);
    bar([refPhase(k,:);rePhase(k,:);rzPhase(k,:)]');
    grid on;
    ylim([-180 180]);
    yticks(-180:60:180);
    set(gca,'XTickLabel',names);
    title(sprintf('%g Hz',freq(k)));

    if k == 1
        ylabel('Phase (degrees)');
    end
end

legend('Direct 0.1 \mus','Mixed Richardson','Rezonance',...
    'Location','best');

sgtitle('Average GFMI admittance phase comparison');

exportgraphics(gcf,fullfile(figureDir,...
    'average_gfmi_rezonance_phase.png'),'Resolution',300);
figure;
bar(freq,[matrixRE matrixRZ]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Full-matrix error (%)');
legend('Mixed Richardson','Rezonance','Location','best');
title('Error against direct 0.1 \mus reference');
set(gca,'YScale','log');
exportgraphics(gcf,fullfile(figureDir,...
    'average_gfmi_rezonance_error.png'),'Resolution',300);

%% Convert Python-style complex values
function z = toComplex(x)
    x = replace(string(x),"j","i");
    z = arrayfun(@(s) str2num(char(s)),x);
end