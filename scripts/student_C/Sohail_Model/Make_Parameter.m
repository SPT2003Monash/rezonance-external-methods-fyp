System=struct();
syms s

%% Base ratings
System.S_base = 2e6;
System.V_LL   = 0.69e3;
System.f_base = 50;
System.w_base = 2*pi*System.f_base;
System.wn = System.w_base;

%% Per-unit base values
System.V_base = (System.V_LL/sqrt(3));
System.I_base = (System.S_base/3)/(System.V_LL/sqrt(3));
System.Z_base = System.V_base/System.I_base;
System.L_base = System.Z_base/System.w_base;
System.C_base = 1/(System.w_base*System.Z_base);

%% Grid model
System.SCR = 5;
System.XR  = 7;
System.rg_pu = 1/(System.SCR*sqrt(1+System.XR^2));
System.lg_pu = System.XR*System.rg_pu;

%% Simulation settings (Not relevant to analytical)
System.f_sample = 40e3;
System.f_step   = 100e3;

%% IBR Parameters
IBR=struct();
IBR.P_pu = 0;
IBR.Q_pu = 0;

%% Most changing parameters %%
%For VI
IBR.lv_pu=0.1;
IBR.rv_pu=0;

%% IBR Physical Parameters
IBR.fsw=5e3;
IBR.td = (0.5/IBR.fsw)+(1/System.f_sample);
IBR.Vdc = 1.45e3;
IBR.w=System.w_base;

%% Filter
IBR.lf1_pu = 0.1;
IBR.rf1_pu = 0.001;
IBR.cf_pu  = 0.05;
IBR.rd_pu  = 0;

%% Conversion using System bases:
IBR.lf1 = IBR.lf1_pu * System.L_base;
IBR.rf1 = IBR.rf1_pu * System.Z_base;
IBR.cf  = IBR.cf_pu  * System.C_base;

%% Current controller
w_cc = 1500;
IBR.k_pi_pu = (IBR.lf1 * w_cc) * (System.I_base/System.V_base);
IBR.k_ii_pu = (IBR.rf1 * w_cc) * (System.I_base/System.V_base);


%% Feed+forwards
IBR.beta_i=0.8;
IBR.beta_v=0.8;

%% Voltage controller
IBR.k_pv_pu=3.3;
IBR.k_iv_pu=176;

%% Reactive power control
IBR.KpQ_pu = 0*0.05;
IBR.KiQ_pu =0*5;

%% APC
IBR.mp= 0.017;% 1% setting
IBR.RoCoF_pu=0.5/System.f_base; %1Hz/s
IBR.H = 1/(2*IBR.RoCoF_pu);
IBR.H = 0.1*IBR.H;
IBR.Dp = 1/IBR.mp;

% save('IBR_Parameters.mat',"IBR");

clear var;