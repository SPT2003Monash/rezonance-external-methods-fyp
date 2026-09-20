% Name: Pahan Kitthangodage
% Student ID: 33597863
% Purpose: Compares switching-GFMI Richardson accuracy and PSCAD runtime using Pareto plots.

clear; clc; close all;

%% DATA 

methods = {'12.5 \mus direct','6.25 \mus direct','Richardson'};
freqs   = [1 5 10];

% Individual PSCAD runs, as listed in the doc: [12.5us d, 12.5us q, 6.25us d, 6.25us q]
runs12p5_d = [60.28 60.40 60.87];
runs12p5_q = [61.35 60.89 61.07];
runs6p25_d = [75.37 78.71 78.17];
runs6p25_q = [79.87 80.07 78.99];

runtime12p5 = runs12p5_d + runs12p5_q;                 % coarse-only, d+q
runtime6p25 = runs6p25_d + runs6p25_q;                 % fine-only, d+q
runtimeRE   = runtime12p5 + runtime6p25;               % Richardson needs both

runtime = [runtime12p5; runtime6p25; runtimeRE]';      % rows = freq, cols = methods

% Full-matrix Frobenius error vs 0.1us reference
% [12.5 us direct, 6.25 us direct, p=2 Richardson]

errorVal = [0.141832 0.039235 0.009786;    % 1 Hz
            2.044435 0.549531 0.201158;    % 5 Hz
            3.826091 1.418587 1.358913];   % 10 Hz
%% PLOT EACH FREQUENCY 

for i = 1:numel(freqs)
    plotPareto(runtime(i,:),errorVal(i,:),methods,freqs(i));
end

%% LOCAL FUNCTION 

function plotPareto(runtime,errorVal,methods,freq)
% A method is Pareto-efficient if no other method has <= runtime AND
% <= error, with a strict improvement in at least one.

    n = numel(runtime);
    isPareto = true(1,n);
    for i = 1:n
        for j = 1:n
            if i==j, continue; end
            if runtime(j)<=runtime(i) && errorVal(j)<=errorVal(i) && ...
               (runtime(j)<runtime(i) || errorVal(j)<errorVal(i))
                isPareto(i) = false;
                break;
            end
        end
    end

    figure; hold on; grid on; box on;

    scatter(runtime(~isPareto),errorVal(~isPareto),70,[0.6 0.6 0.6],'filled');
    scatter(runtime(isPareto),errorVal(isPareto),90,[0.85 0.10 0.10],'filled');

    [prRuntime,order] = sort(runtime(isPareto));
    prError = errorVal(isPareto); prError = prError(order);
    plot(prRuntime,prError,'-','Color',[0.85 0.10 0.10],'LineWidth',1.5);

    xlim([min(runtime)-5, max(runtime)+15]);
    ylim([min(errorVal)*0.85, max(errorVal)*1.35]);
    set(gca,'YScale','log');

    [~,orderByX] = sort(runtime);
    for idx = 1:n
        k = orderByX(idx);
        closeToPrev = idx>1 && abs(runtime(k)-runtime(orderByX(idx-1))) < (0.02*(max(runtime)-min(runtime)));
        if closeToPrev
            vAlign = 'top'; dy = 0.94;
        else
            vAlign = 'bottom'; dy = 1.06;
        end
        text(runtime(k)+2,errorVal(k)*dy,methods{k},'FontSize',9,'VerticalAlignment',vAlign);
    end

    xlabel('PSCAD serial runtime, d + q runs (s)');
    ylabel('Full-matrix error vs 0.1 \mus reference (%)');
    title(sprintf('%g Hz — Switching GFMI Accuracy vs Runtime',freq));
    legend({'Dominated','Pareto-efficient','Pareto front'},'Location','best');
    hold off;

end
