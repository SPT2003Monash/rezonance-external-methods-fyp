% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Average-GFMI accuracy/runtime Pareto plots

clear; clc; close all;

methods = {'50 \mus direct','25 \mus direct','12.5 \mus direct', ...
           'p=1 RE','p=2 RE','Mixed RE'};
freqs = [1 5 10];

% Total PSCAD runtime = d-perturbation + q-perturbation
runtime = [10.83 18.77 32.86 29.60 29.60 62.46;   % 1 Hz
           10.54 18.50 33.30 29.04 29.04 62.34;   % 5 Hz
           10.41 18.44 33.56 28.85 28.85 62.41];  % 10 Hz

% Full-matrix error relative to direct 0.1 us reference
errorVal = [0.315689197 0.157496035 0.078295238 ...
            0.001564820 0.104747888 0.001188091;   % 1 Hz
            1.570754440 0.782189658 0.389144758 ...
            0.020382463 0.519185114 0.003077358;   % 5 Hz
            2.619480691 1.281278190 0.632064201 ...
            0.072444472 0.836859052 0.004289492];  % 10 Hz

% 0.1 us d + q runtimes
refRuntime = [4510 4563 4588];

matlabProcessingTime = mean([0.395127 0.413888 0.374571]);

fprintf('\n=== UPDATED AVERAGE GFMI PARETO ANALYSIS ===\n');
fprintf('Average MATLAB processing time = %.3f s\n', ...
        matlabProcessingTime);

for i = 1:numel(freqs)
    fprintf('\n--- %g Hz ---\n',freqs(i));

    for k = 1:numel(methods)
        fprintf('%-18s runtime = %6.2f s   error = %.9f %%\n', ...
            methods{k},runtime(i,k),errorVal(i,k));
    end

    fprintf('0.1 us reference runtime = %.2f s\n',refRuntime(i));
end

fprintf('\n=== SPEED-UP RELATIVE TO 0.1 us REFERENCE ===\n');

for i = 1:numel(freqs)
    fprintf('\n%g Hz\n',freqs(i));

    for k = 1:numel(methods)
        fprintf('%-18s = %.2fx faster\n',methods{k}, ...
            refRuntime(i)/runtime(i,k));
    end
end

fprintf('\n=== BEST RICHARDSON METHOD ===\n');

for i = 1:numel(freqs)
    [bestErr,j] = min(errorVal(i,4:6));
    k = j + 3;

    fprintf('%g Hz: %s | error = %.9f %% | runtime = %.2f s\n', ...
        freqs(i),methods{k},bestErr,runtime(i,k));
end

for i = 1:numel(freqs)
    plotPareto(runtime(i,:),errorVal(i,:),methods, ...
               refRuntime(i),freqs(i));
end

function plotPareto(runtime,errorVal,methods,referenceRuntime,freq)

n = numel(runtime);
isPareto = true(1,n);

for i = 1:n
    for j = 1:n
        if i ~= j && runtime(j)<=runtime(i) && ...
           errorVal(j)<=errorVal(i) && ...
           (runtime(j)<runtime(i) || errorVal(j)<errorVal(i))

            isPareto(i) = false;
            break
        end
    end
end

figure; hold on; grid on; box on;

scatter(runtime(~isPareto),errorVal(~isPareto),70, ...
        [0.6 0.6 0.6],'filled');

scatter(runtime(isPareto),errorVal(isPareto),90, ...
        [0.85 0.10 0.10],'filled');

[xp,order] = sort(runtime(isPareto));
yp = errorVal(isPareto);
plot(xp,yp(order),'-','Color',[0.85 0.10 0.10], ...
     'LineWidth',1.5);

xline(referenceRuntime,'--', ...
    sprintf('0.1 \\mus reference (%.0f s)',referenceRuntime));

set(gca,'XScale','log','YScale','log');

for k = 1:n
    text(runtime(k)*1.04,errorVal(k)*1.08,methods{k}, ...
         'FontSize',9);
end

xlabel('PSCAD serial runtime, d + q runs (s)');
ylabel('Full-matrix error vs 0.1 \mus reference (%)');
title(sprintf('OP1, %g Hz — Accuracy vs Runtime',freq));

legend({'Dominated','Pareto-efficient','Pareto front'}, ...
       'Location','best');

hold off
end