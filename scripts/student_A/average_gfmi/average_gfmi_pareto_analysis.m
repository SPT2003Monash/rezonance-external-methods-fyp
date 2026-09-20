% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares average-GFMI accuracy and PSCAD runtime using Pareto plots.

clear; clc; close all;

%% DATA -------------------------------------------------------------------
% Runtime columns: 50 direct | 25 direct | 12.5 direct | p=1 RE (50+25) |
%                  p=2 RE (50+25) | Mixed RE (50+25+12.5)
% Each direct runtime = d-perturbation + q-perturbation PSCAD run time.

methods = {'50 \mus direct','25 \mus direct','12.5 \mus direct','p=1 RE','p=2 RE','Mixed RE'};
freqs   = [1 5 10];

runtime = [ 5.31  9.61 16.51 14.92 14.92 31.43;   % 1 Hz
            5.68  9.55 16.77 15.23 15.23 32.00;   % 5 Hz
            5.38  9.44 16.63 14.82 14.82 31.45];  % 10 Hz

errorVal = [0.301121053 0.145201983 0.067278550 0.012588306 0.093259350 0.012635929;   % 1 Hz
            1.502475921 0.729145561 0.343558495 0.048325252 0.471270592 0.041317892;   % 5 Hz
            2.420598896 1.156714664 0.541596658 0.112272025 0.737177889 0.062654918];  % 10 Hz

refRuntime = [95.60 96.38 96.54];   % 2 us reference runtime (zero error by definition, not plotted)
matlabProcessingTime = 0.521837;    % average of final 1, 5 and 10 Hz MATLAB runs

fprintf('\n=== VALIDATED AVERAGE GFMI PARETO ANALYSIS ===\n');
fprintf('Average MATLAB processing time = %.3f s\n', matlabProcessingTime);

%% RUNTIME / ERROR TABLE ----------------------------------------------------

for i = 1:numel(freqs)
    fprintf('\n--- %g Hz ---\n', freqs(i));
    for k = 1:numel(methods)
        fprintf('%-18s runtime = %6.2f s   error = %.9f %%\n', methods{k}, runtime(i,k), errorVal(i,k));
    end
    fprintf('2 us reference  runtime = %.2f s\n', refRuntime(i));
end

%% SPEED-UP RELATIVE TO 2 us REFERENCE ------------------------------------

fprintf('\n=== SPEED-UP RELATIVE TO DIRECT 2 us REFERENCE ===\n');
for i = 1:numel(freqs)
    fprintf('\n%g Hz\n', freqs(i));
    for k = 1:numel(methods)
        fprintf('%-18s = %.3fx faster\n', methods{k}, refRuntime(i)/runtime(i,k));
    end
end

%% BEST RICHARDSON METHOD PER FREQUENCY --------------------------------------

fprintf('\n=== BEST RICHARDSON METHOD AT EACH FREQUENCY ===\n');
REidx = [4 5 6];
for i = 1:numel(freqs)
    [bestErr,idxBest] = min(errorVal(i,REidx));
    bestMethodIndex = REidx(idxBest);
    fprintf('%g Hz: %s | error = %.9f %% | runtime = %.2f s\n', ...
        freqs(i), methods{bestMethodIndex}, bestErr, runtime(i,bestMethodIndex));
end

%% PLOT EACH FREQUENCY -------------------------------------------------------

for i = 1:numel(freqs)
    plotPareto(runtime(i,:), errorVal(i,:), methods, refRuntime(i), freqs(i));
end

%% ============================ LOCAL FUNCTION ================================

function plotPareto(runtime, errorVal, methods, referenceRuntime, freq)
    % A method is Pareto-efficient if no other method has <= runtime AND
    % <= error, with a strict improvement in at least one.
    n = numel(runtime);
    isPareto = true(1,n);
    for i = 1:n
        for j = 1:n
            if i == j; continue; end
            if runtime(j)<=runtime(i) && errorVal(j)<=errorVal(i) && ...
               (runtime(j)<runtime(i) || errorVal(j)<errorVal(i))
                isPareto(i) = false;
                break;
            end
        end
    end

    figure; hold on; grid on; box on;

    % Dominated methods: grey. Pareto-efficient: highlighted + connected.
    scatter(runtime(~isPareto), errorVal(~isPareto), 70, [0.6 0.6 0.6], 'filled');
    scatter(runtime(isPareto),  errorVal(isPareto),  90, [0.85 0.10 0.10], 'filled');

    [prRuntime,order] = sort(runtime(isPareto));
    prError = errorVal(isPareto); prError = prError(order);
    plot(prRuntime, prError, '-', 'Color',[0.85 0.10 0.10], 'LineWidth',1.5);

    xline(referenceRuntime, '--', sprintf('2 \\mus reference (%.1f s)', referenceRuntime), ...
        'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal');

    xlim([min(runtime)-1, max([runtime referenceRuntime])+3]);
    ylim([min(errorVal)*0.85, max(errorVal)*1.35]);
    set(gca,'YScale','log');

    % Offset each method's text label; nudge below the point if two
    % methods share nearly the same runtime, to avoid overlapping labels.
    [~,orderByX] = sort(runtime);
    for idx = 1:n
        k = orderByX(idx);
        closeToPrev = idx>1 && abs(runtime(k)-runtime(orderByX(idx-1))) < 0.5;
        if closeToPrev; vAlign='top'; dy=0.94; else; vAlign='bottom'; dy=1.06; end
        text(runtime(k)+0.3, errorVal(k)*dy, methods{k}, 'FontSize',9, 'VerticalAlignment',vAlign);
    end

    xlabel('PSCAD serial runtime, d + q runs (s)');
    ylabel('Full-matrix error vs 2 \mus reference (%)');
    title(sprintf('%g Hz — Accuracy vs Runtime Pareto Plot', freq));
    legend({'Dominated','Pareto-efficient','Pareto front'}, 'Location','best');
    hold off;
end