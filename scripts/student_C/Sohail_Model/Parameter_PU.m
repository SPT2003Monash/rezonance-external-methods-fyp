%% Base ratings
S_base = 2e6;
V_LL   = 0.69e3;
f_base = 50;
w_base = 2*pi*f_base;
wn = w_base;

%% Per-unit base values
V_base = (V_LL/sqrt(3));
I_base = (S_base/3)/(V_LL/sqrt(3));
Z_base = V_base/I_base;
L_base = Z_base/w_base;
C_base = 1/(w_base*Z_base);

%% Grid model
SCR = 5;
XR  = 7;
rg_pu = 1/(SCR*sqrt(1+XR^2));
lg_pu = XR*rg_pu;

%% Simulation settings (Not relevant to analytical)
f_sample = 40e3;
f_step   = 100e3;

%% IBR Parameters
P_pu = 0;
Q_pu = 0;

%% Most changing parameters %%
%For VI
lv_pu = 0.1;
rv_pu = 0;

%% IBR Physical Parameters
fsw = 5e3;
td  = (0.5/fsw)+(1/f_sample);
Vdc = 1.45e3;
w   = w_base;

%% Filter
lf1_pu = 0.1;
rf1_pu = 0.002;
cf_pu  = 0.04;
rd_pu  = 0.01;

%% Conversion using bases:
lf1 = lf1_pu * L_base;
rf1 = rf1_pu * Z_base;
cf  = cf_pu  * C_base;
rd = rd_pu * Z_base;
%% Current controller
w_cc = 1500;
k_pi_pu = 0.318;%(lf1 * w_cc) * (I_base/V_base);
k_ii_pu = 2;%(rf1 * w_cc) * (I_base/V_base);

%% Feed-forwards
beta_i = 0.837;
beta_v = 0.5;

%% Voltage controller
k_pv_pu = 0.76;
k_iv_pu = 292;

%% Reactive power control
KpQ_pu = 0.05;
KiQ_pu = 5;

%% APC
mp = 0.017; % 1% setting
RoCoF_pu = 0.5/f_base; %1Hz/s
H = 1/(2*RoCoF_pu);
H = 0.1*H;
Dp = 1/mp;

save('Parameters_Sohail.mat');