%% Power Injection Setpoints
P_inj = 1.6;                 % Active power injection (pu)
Q_inj = 0.2;                 % Reactive power injection (pu)

%% General System Parameters
f_base = 50;                 % Base frequency (Hz)
wn = 2 * pi * f_base;        % Angular base frequency (rad/s)
w = wn;
% f_sample = 20e3;             % Sample frequency 20kHz
% f_step = 100e3;              % Step frequency 100kHz

%% Inverter Ratings
Vdc = 1.2e3;                 % DC-link voltage (V)
V_LL = 0.69e3;               % Line-to-line RMS voltage (V)
S_base = 2e6;                % Rated power of IBR (1 MW)
Zb = V_LL^2 / S_base;        % Base impedance (Ohms)
Lb = Zb / wn;                % Base inductance (H)
Cb = 1/(Zb * wn);


%% Grid Model Parameters
SCR = 5;                     % Short Circuit Ratio (SCR)
XR = 5;                      % X/R Ratio
rgpu = 1 / (SCR * sqrt(1 + XR^2)); % Grid resistance in pu
lgpu = XR * rgpu;            % Grid inductance in pu
rg = rgpu * Zb;              % Grid resistance (Ohms)
lg = lgpu * Lb;              % Grid inductance (H)

%% Switching and Control Delays
fsw = 5e3;                   % Switching frequency (Hz)
td = 1.5 / fsw;              % Delay accounting for inner control loops (s)

%% LC Filter Parameters (TO change) 
lf1 = 0.1 * Lb;                % Inverter-side inductor (H)
rf1 = 0.002 * Zb;             % Inverter-side resistance (Ohms)


cf = 0.04 * Cb;                 % Filter capacitor (F)
rd = 0.01 * Zb;                  % Damping resistor (Ohms)

L_t = lf1 + lg;
f_res = (1/(2*pi)) * sqrt((lf1 + lg) / (lf1 *(lg)* cf));
f_ares = (1/(2*pi)) * sqrt(1 / (lf1 * cf));

%% Current Control Parameters (To change and modify)
Tic = 2 * td;                % Current control time constant (s)
k_pi = 0.318*Zb;%lf1 / Tic;            % Proportional gain of current controller
k_ii = 2*Zb;%rf1 / Tic;            % Integral gain of current controller
beta_v = 0.5;                % Voltage feedforward gain factor

%% Voltage Control Parameters 
wm = 2000;                   % Overriding default setting (rad/s)
k_pv = 0.764/Zb;%cf * wm;              % Proportional gain of voltage controller
k_iv = 292/Zb;%Tic * cf * wm^3;      % Integral gain of voltage controller
beta_i = 0.837;                % Current feedforward gain factor

%% Active Power Control (APC) - Droop + LPF 
mp = 0.5/(S_base); %0.5 scaled to real values. 
wc_p = 200;              % Power LPF cutoff frequency (~5 Hz)

%% Reactive Power Control (RPC) - Droop + LPF
% In amplitude-invariant dq transform, nominal Vd is the peak phase voltage.
V_set = sqrt(2/3) * V_LL;    % Nominal dq-frame voltage (approx 563.38 V)
nq = 0.05*V_set/S_base;%V_drop / S_base;        % Reactive power droop gain (Volts / VAr)
wc_q = 100;                 % Reactive power LPF cutoff frequency (~5 Hz)

%% Virtual Impedance 
% Setting virtual inductance to approx 10% of base impedance to ensure P/Q decoupling
lv = 0.03 * Lb;               % Virtual Inductance (H)
rv = 0.003*Zb;                      % Virtual Resistance (Ohms)

R_vir = rv;                  % Passed to State-Space Model
X_vir = wn * lv;             % Converted to Virtual Reactance for State-Space Model

%% Save Parameters for Simulink
save('Parameters_LC.mat', 'f_base', 'f_sample','f_step', 'wn', 'w','V_LL', 'S_base', 'Zb', 'Vdc', 'Lb', ...
     'fsw','f_res', 'f_ares', 'td', 'lf1', 'rf1','cf', 'rd', 'SCR', 'XR', 'rgpu', 'lgpu', 'rg', 'lg', 'k_pi', 'k_ii', ...
     'beta_v', 'k_pv','k_iv','beta_i','Tic', 'P_inj', 'Q_inj', ...
     'mp', 'wc_p', 'V_set', 'nq', 'wc_q', 'R_vir', 'X_vir');
