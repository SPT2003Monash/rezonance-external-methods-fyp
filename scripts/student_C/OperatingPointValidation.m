%% Corrected Power Flow & d-q Transformation Script
mpc = SMIB_PowerFlow(rg, lg);  
PowerFlow_results = runpf(mpc);

% Extract bus 1 angle (in degrees) and voltage magnitude (p.u)
Angle = PowerFlow_results.bus(1, 9); % Load angle in degrees (e.g. +12 deg)
Voltage_mag = PowerFlow_results.bus(1, 8);    % Voltage magnitude in p.u.

% 1. Calculate Line-to-Neutral Voltage Phasors in the Grid Reference Frame (Vg = 0 deg)
Angle_rad = deg2rad(Angle); % Positive leading angle
V_LN_RMS = (V_LL / sqrt(3));

Vc_phasor = Voltage_mag * V_LN_RMS * exp(1j * Angle_rad); % Vc leads Vg
Vg_phasor = V_LN_RMS * exp(1j * 0);                       % Grid aligned at 0 deg

% 2. Calculate Complex Branch Currents in Grid Frame (RMS)
Zt = (rf2 + rg) + 1j * 2 * pi * 50 * (lf2 + lg);
Z2 = rf2 + 1j * 2 * pi * 50 * lf2;

I2_phasor = (Vc_phasor - Vg_phasor) / Zt; % Grid-side current phasor (RMS)
Vpcc_phasor = Vc_phasor - I2_phasor * Z2;

% 3. Transform Phasors to Inverter d-q Frame (Aligned with Vc)
% To align frame with Vc, multiply by exp(-1j * Angle_rad) and sqrt(2) for peak
Vc_dq0 = Vc_phasor * exp(-1j * (1.75-pi/2)) * sqrt(2); % Corrected minus sign!
Vc_d0  = real(Vc_dq0); % Should equal Voltage_mag * V_LN_RMS * sqrt(2)
Vc_q0  = imag(Vc_dq0); % Should equal exactly 0.0

I2_dq0 = I2_phasor * exp(-1j * (1.75-pi/2)) * sqrt(2);
i2_d0  = real(I2_dq0) % Will now output ~1859 A!
i2_q0  = imag(I2_dq0)

% Store operating points
delta0  = Angle_rad;    
Vpcc_D0 = real(Vpcc_phasor * exp(-1j * Angle_rad) * sqrt(2));
Vpcc_Q0 = imag(Vpcc_phasor * exp(-1j * Angle_rad) * sqrt(2));
%%
%mathematically correct here. 1.6MW and 0.2MVAR
%operating point validation. 
realpower = 3/2*(Vc_d0*i2_d0 + Vc_q0*i2_q0)
reactive = -3/2*(Vc_q0*i2_d0-Vc_d0*i2_q0)
%%
% 4. Calculate Missing Steady-State Phasors (I1 and Vcf)
omega = 2 * pi * 50; % Nominal angular frequency

% Capacitor current (Ic) using Vc_phasor and filter capacitance (cf)
Ic_phasor = Vc_phasor * (1j * omega * cf);

% Inverter-side current (I1) using KCL: I1 = I2 + Ic
I1_phasor = I2_phasor + Ic_phasor;

% Inverter terminal voltage (Vcf) using KVL: Vcf = Vc + I1 * Z1
Z1 = rf1 + 1j * omega * lf1;
Vcf_phasor = Vc_phasor + I1_phasor * Z1;

% 5. Transform Missing Phasors to Inverter d-q Frame
% Transform I1
I1_dq0 = I1_phasor * exp(-1j * PSCAD_angle) * sqrt(2);
i1_d0  = real(I1_dq0);
i1_q0  = imag(I1_dq0);

% Transform Vcf
Vcf_dq0 = Vcf_phasor * exp(-1j * PSCAD_angle) * sqrt(2);
Vcf_d0  = real(Vcf_dq0);
Vcf_q0  = imag(Vcf_dq0);

%Simulatino values. 
i2dq = [1859 125]
% %%
% % 1. Define frequencies and remaining impedances
% omega = 2 * pi * 50;
% Z_f1 = rf1 + 1j * omega * lf1;  % Inverter-side filter impedance
% Z_c  = rd + 1 / (1j * omega * cf); % Capacitor branch impedance
% 
% % 2. Calculate Capacitor Current (KCL)
% Ic_phasor = Vc_phasor / Z_c;
% 
% % 3. Calculate Inverter Current I1 (This fixes your 4642 A error!)
% I1_phasor = I2_phasor + Ic_phasor;
% 
% % 4. Calculate Inverter Internal Voltage E_inv (KVL)
% E_inv_phasor = Vc_phasor + (I1_phasor * Z_f1);
% 
% % 5. Extract the Angles (in radians)
% theta_inv = angle(E_inv_phasor); % Absolute angle of the internal voltage
% theta_Vc = angle(Vc_phasor);     % Absolute angle of the capacitor voltage
% 
% % 6. Calculate the control load angle delta
% delta = theta_inv - theta_Vc;
% 
% %%
% % 1. Convert RMS phasor to Peak magnitude
% I1_peak = I1_phasor * sqrt(2);
% 
% % 2. Rotate from the global XY stationary frame into the Inverter's rotating d-q frame
% % (Multiplying by exp(-1j * theta_inv) aligns the d-axis with E_inv)
% I1_dq = I1_peak * exp(-1j * theta_inv);
% 
% % 3. Extract the matched d and q components (Validated within PSCAD). 
% i1_d0 = real(I1_dq)
% i1_q0 = -imag(I1_dq)
% 
% %% Match I2 to PSCAD
% % 1. Convert RMS phasor to Peak magnitude
% I2_peak = I2_phasor * sqrt(2);
% 
% % 2. Rotate into the Inverter's rotating d-q frame (aligned with E_inv)
% I2_dq_inverter_frame = I2_peak * exp(-1j*theta_inv);
% 
% % 3. Extract components and flip q-axis to match PSCAD lagging convention
% i2_d0 = 1859
% i2_q0 = 119


