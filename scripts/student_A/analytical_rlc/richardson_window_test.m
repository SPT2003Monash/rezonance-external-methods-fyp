% Name : Pahan Kitthangodage
% Student ID : 33597863
% Purpose: Tests Richardson refinement using synthetic short-, medium- and long-window dq impedance estimates.

clear; clc; close all;

%% 1. Frequency vector
f = logspace(-1, 2, 200);   % 0.1 Hz to 100 Hz
w = 2*pi*f;
s = 1j*w;

%% 2. Parameters for synthetic impedance model
R = 0.5;        % Resistance
L = 0.02;       % Inductance

% Create synthetic 2x2 dq impedance matrix at each frequency
Z_true = zeros(2,2,length(f));

for k = 1:length(f)
    Z_base = R + s(k)*L;

    % Synthetic 2x2 dq impedance matrix used for method verification
    % Diagonal terms represent main d-axis and q-axis impedance behaviour.
    % Off-diagonal terms represent small dq coupling.
    Z_true(:,:,k) = [ Z_base,        0.05*Z_base;
                    -0.03*Z_base,   1.2*Z_base ];
end

%% 3. Week 10 scan-window settings
% Larger h means lower quality / shorter measurement window / larger error.
% Smaller h means higher quality / longer measurement window / smaller error.

h_short  = 0.4;    % Short-window scan: fastest but least accurate
h_medium = 0.2;    % Medium-window scan: moderate accuracy
h_long   = 0.1;    % Long-window scan: slowest but best accuracy

p = 2;             % Assumed convergence order for Richardson formula

% Refinement ratios
r_SM = h_short / h_medium;   % short to medium
r_ML = h_medium / h_long;    % medium to long
r_SL = h_short / h_long;     % short to long

%% 4. Noise settings
rng(1);   % Makes random noise repeatable

% Different noise levels imitate different scan-window qualities.
% Short window has more noise, long window has less noise.
noise_short  = 0.003;
noise_medium = 0.0015;
noise_long   = 0.0007;

%% 5. Preallocate matrices
Z_short  = zeros(size(Z_true));
Z_medium = zeros(size(Z_true));
Z_long   = zeros(size(Z_true));

Z_R_SM = zeros(size(Z_true));   % Richardson using short + medium
Z_R_ML = zeros(size(Z_true));   % Richardson using medium + long
Z_R_SL = zeros(size(Z_true));   % Richardson using short + long

% Error arrays
err_short  = zeros(1,length(f));
err_medium = zeros(1,length(f));
err_long   = zeros(1,length(f));

err_R_SM = zeros(1,length(f));
err_R_ML = zeros(1,length(f));
err_R_SL = zeros(1,length(f));

% Acceptance flags for Richardson estimates
accept_R_SM = false(1,length(f));
accept_R_ML = false(1,length(f));
accept_R_SL = false(1,length(f));

%% 6. Generate synthetic scan estimates and apply Richardson refinement
for k = 1:length(f)

    % Artificial complex error coefficient
    % This represents the systematic scan error that reduces as h becomes smaller.
    C = 0.15 * Z_true(:,:,k) * exp(1j*pi/6);

    % Random complex noise for each scan-window case
    noiseS = noise_short * norm(Z_true(:,:,k),'fro') * ...
        (randn(2,2) + 1j*randn(2,2));

    noiseM = noise_medium * norm(Z_true(:,:,k),'fro') * ...
        (randn(2,2) + 1j*randn(2,2));

    noiseL = noise_long * norm(Z_true(:,:,k),'fro') * ...
        (randn(2,2) + 1j*randn(2,2));

    % Synthetic scan estimates
    Z_short(:,:,k)  = Z_true(:,:,k) + C*h_short^p  + noiseS;
    Z_medium(:,:,k) = Z_true(:,:,k) + C*h_medium^p + noiseM;
    Z_long(:,:,k)   = Z_true(:,:,k) + C*h_long^p   + noiseL;

    % Richardson refinement estimates
    Z_R_SM(:,:,k) = (r_SM^p * Z_medium(:,:,k) - Z_short(:,:,k)) / (r_SM^p - 1);
    Z_R_ML(:,:,k) = (r_ML^p * Z_long(:,:,k)   - Z_medium(:,:,k)) / (r_ML^p - 1);
    Z_R_SL(:,:,k) = (r_SL^p * Z_long(:,:,k)   - Z_short(:,:,k)) / (r_SL^p - 1);

    % Relative matrix error using Frobenius norm
    err_short(k) = norm(Z_short(:,:,k)  - Z_true(:,:,k),'fro') / norm(Z_true(:,:,k),'fro');
    err_medium(k)= norm(Z_medium(:,:,k) - Z_true(:,:,k),'fro') / norm(Z_true(:,:,k),'fro');
    err_long(k)  = norm(Z_long(:,:,k)   - Z_true(:,:,k),'fro') / norm(Z_true(:,:,k),'fro');

    err_R_SM(k) = norm(Z_R_SM(:,:,k) - Z_true(:,:,k),'fro') / norm(Z_true(:,:,k),'fro');
    err_R_ML(k) = norm(Z_R_ML(:,:,k) - Z_true(:,:,k),'fro') / norm(Z_true(:,:,k),'fro');
    err_R_SL(k) = norm(Z_R_SL(:,:,k) - Z_true(:,:,k),'fro') / norm(Z_true(:,:,k),'fro');

    %% 7. Simple Richardson acceptance / rejection rule
    % Since this is synthetic data, we can compare to Z_true directly.
    % In real PSCAD/Rezonance data, Z_true will not be available.
    %
    % Accept Richardson only if it improves compared with the better input scan.

    accept_R_SM(k) = err_R_SM(k) < err_medium(k);
    accept_R_ML(k) = err_R_ML(k) < err_long(k);
    accept_R_SL(k) = err_R_SL(k) < err_long(k);

end

%% 8. Plot relative error comparison
figure('Units','centimeters','Position',[2 2 8.8 6]);

loglog(f, err_short,  'LineWidth', 1.1); hold on;
loglog(f, err_medium, 'LineWidth', 1.1);
loglog(f, err_long,   'LineWidth', 1.1);

loglog(f, err_R_SM, '--', 'LineWidth', 1.1);
loglog(f, err_R_ML, '--', 'LineWidth', 1.1);
loglog(f, err_R_SL, '--', 'LineWidth', 1.1);

grid on;
box on;

xlabel('Frequency (Hz)', 'FontSize', 8);
ylabel('Relative Error', 'FontSize', 8);

title('Richardson Refinement', ...
      'FontSize', 8, 'FontWeight', 'bold');

legend('Short', 'Medium', 'Long', ...
       'R: S+M', 'R: M+L', 'R: S+L', ...
       'Location', 'southwest', ...
       'FontSize', 5.5);

set(gca, 'FontSize', 8);
set(gcf, 'Color', 'w');

exportgraphics(gcf, 'fig1_richardson_scan_windows.png', 'Resolution', 300);
%% 9. Plot accepted/rejected Richardson points
figure;
semilogx(f, accept_R_SM, 'LineWidth', 1.5); hold on;
semilogx(f, accept_R_ML, 'LineWidth', 1.5);
semilogx(f, accept_R_SL, 'LineWidth', 1.5);

grid on;
xlabel('Frequency (Hz)');
ylabel('Accepted = 1, Rejected = 0');
legend('Richardson short + medium', ...
       'Richardson medium + long', ...
       'Richardson short + long', ...
       'Location', 'best');

title('Week 10: Simple Richardson Acceptance/Rejection Check');
ylim([-0.1 1.1]);

%% 10. Print summary table
fprintf('\n================ Week 10 Richardson Window Test Summary ================\n');

fprintf('\nMedian Relative Error:\n');
fprintf('Short-window scan:              %.4e\n', median(err_short));
fprintf('Medium-window scan:             %.4e\n', median(err_medium));
fprintf('Long-window scan:               %.4e\n', median(err_long));
fprintf('Richardson short + medium:      %.4e\n', median(err_R_SM));
fprintf('Richardson medium + long:       %.4e\n', median(err_R_ML));
fprintf('Richardson short + long:        %.4e\n', median(err_R_SL));

fprintf('\nMaximum Relative Error:\n');
fprintf('Short-window scan:              %.4e\n', max(err_short));
fprintf('Medium-window scan:             %.4e\n', max(err_medium));
fprintf('Long-window scan:               %.4e\n', max(err_long));
fprintf('Richardson short + medium:      %.4e\n', max(err_R_SM));
fprintf('Richardson medium + long:       %.4e\n', max(err_R_ML));
fprintf('Richardson short + long:        %.4e\n', max(err_R_SL));

%% 11. Print acceptance percentage
accept_percent_SM = 100 * sum(accept_R_SM) / length(f);
accept_percent_ML = 100 * sum(accept_R_ML) / length(f);
accept_percent_SL = 100 * sum(accept_R_SL) / length(f);

fprintf('\nRichardson Acceptance Percentage:\n');
fprintf('Short + medium accepted:        %.1f %% of frequency points\n', accept_percent_SM);
fprintf('Medium + long accepted:         %.1f %% of frequency points\n', accept_percent_ML);
fprintf('Short + long accepted:          %.1f %% of frequency points\n', accept_percent_SL);

%% 12. Simple conclusion printed by code
fprintf('\nConclusion:\n');

if median(err_R_ML) < median(err_long)
    fprintf('Richardson medium + long improved the median error compared with the long-window scan.\n');
else
    fprintf('Richardson medium + long did NOT improve the median error compared with the long-window scan.\n');
end

if median(err_R_SM) < median(err_medium)
    fprintf('Richardson short + medium improved the median error compared with the medium-window scan.\n');
else
    fprintf('Richardson short + medium did NOT improve the median error compared with the medium-window scan.\n');
end

fprintf('This shows that Richardson refinement should be checked before use, because noise can make the refined estimate worse.\n');
fprintf('========================================================================\n');