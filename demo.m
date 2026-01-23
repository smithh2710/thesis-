clear; clc; close all;


components = {'N2', 'CO2', 'C1', 'C2', 'C3', 'iC4', 'nC4', 'iC5', 'nC5', ...
              'C6', 'C7', 'C8', 'C9', 'C10-C11', 'C12-C13', 'C14-C16', ...
              'C17-C18', 'C19-C21', 'C22-C24', 'C25-C29', 'C30-C37', 'C38-C80'};

% Reference Composition at 175 m (mol%) - Table 7
comp_ref = [0.42; 0.69; 50.04; 7.85; 6.77; 1.04; 3.20; 1.16; 1.55; 1.88; ...
            3.50; 3.75; 2.28; 3.26; 2.59; 2.93; 1.46; 1.65; 1.17; 1.24; ...
            0.96; 0.63];
comp_ref = comp_ref / sum(comp_ref);

% Molecular Weights (g/mol)
M_gmol = [28.014 ; 44.010 ; 16.043 ; 30.070 ; 44.097 ; 58.124 ; 58.124 ; 72.151 ; 72.151 ; 86.178 ; ...
          96.0; 107.0; 121.0; 140.1; 167.6; 204.7; 243.6; 275.3; 317.0; ...
          370.4; 456.8; 640.9];

% Critical Temperatures (K)
Tc = [-147.0; 31.1; -82.5; 32.2; 96.6; 134.9; 152.1; 187.2; 196.4; ...
      234.2; 269.9; 290.9; 315.3; 345.9; 385.0; 432.7; 476.8; 511.1; ...
      553.3; 604.9; 683.3; 847.8] + 273.15;

% Critical Pressures (Pa)
Pc = [33.94; 73.76; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; ...
      29.69; 29.60; 27.86; 25.83; 23.76; 21.58; 19.59; 18.24; 17.54; ...
      16.84; 16.23; 15.62; 14.72] * 1e5;

% Acentric Factors
acentric = [0.040; 0.225; 0.008; 0.098; 0.152; 0.176; 0.193; 0.227; 0.251; ...
            0.296; 0.338; 0.374; 0.420; 0.483; 0.570; 0.685; 0.795; 0.881; ...
            0.984; 1.098; 1.229; 1.159];

% Jhaveri-Youngren Volume Corrections (cm³/mol) - Table 7
c_JY = [-4.23; -1.91; -5.20; -5.79; -6.35; -7.18; -6.49; -6.20; -5.12; ...
        1.42; 8.45; 9.20; 10.32; 11.12; 10.48; 6.57; -1.30; -9.89; ...
        -23.43; -43.30; -79.96; -151.44];

% Critical Volumes (cm³/mol)
R = 8.314472;
Zc_est = 0.2905 - 0.085 * acentric;
Vc = Zc_est .* R .* Tc ./ Pc * 1e6;
Vc_lit = [89.8; 94.0; 99.0; 148.0; 203.0; 263.0; 255.0; 306.0; 304.0; 370.0];
Vc(1:10) = Vc_lit;


n = length(comp_ref);
BIP = zeros(n, n);

% N2 BIPs
BIP(1,2) = -0.0170; BIP(2,1) = -0.0170;
BIP(1,3) = 0.0311;  BIP(3,1) = 0.0311;
BIP(1,4) = 0.0515;  BIP(4,1) = 0.0515;
BIP(1,5) = 0.0852;  BIP(5,1) = 0.0852;
BIP(1,6) = 0.1033;  BIP(6,1) = 0.1033;
BIP(1,7) = 0.0800;  BIP(7,1) = 0.0800;
BIP(1,8) = 0.0922;  BIP(8,1) = 0.0922;
BIP(1,9) = 0.1000;  BIP(9,1) = 0.1000;
for i = 10:n
    BIP(1,i) = 0.0800; BIP(i,1) = 0.0800;
end

% CO2 BIPs
for i = 3:10
    BIP(2,i) = 0.1200; BIP(i,2) = 0.1200;
end
for i = 11:n
    BIP(2,i) = 0.0800; BIP(i,2) = 0.0800;
end


h_ref = 175;
press_ref = 284e5;
temp_ref = 93 + 273.15;

Cp_coeffs = [
    % Component    C1        C2        C3          C4
    31.15,     -0.014,    2.68e-5,   -1.17e-8;    % N2 
    19.79,      0.073,   -5.60e-5,    1.72e-8;    % CO2
    19.25,      0.052,    1.20e-5,   -1.13e-8;    % C1
     5.41,      0.178,   -6.94e-5,    8.71e-9;    % C2
    -4.22,      0.306,   -1.59e-4,    3.21e-8;    % C3
    -1.39,      0.385,   -1.83e-4,    2.90e-8;    % iC4
     9.49,      0.331,   -1.11e-4,   -2.82e-9;    % nC4
    -9.52,      0.507,   -2.73e-4,    5.72e-8;    % iC5
    -3.63,      0.487,   -2.58e-4,    5.30e-8;    % nC5
    -4.41,      0.582,   -3.12e-4,    6.49e-8;    % C6
    -5.15,      0.676,   -3.65e-4,    7.66e-8;    % nC7
    -6.10,      0.771,   -4.20e-4,    8.85e-8;    % nC8
     3.14,      0.677,   -1.93e-4,   -2.98e-8;    % nC9
    25.20,      0.830,   -3.23e-4,    4.06e-8;    % C10-C11
    30.14,      0.993,   -3.87e-4,    4.86e-8;    % C12-C13
    36.83,      1.213,   -4.72e-4,    5.93e-8;    % C14-C16
    43.82,      1.443,   -5.62e-4,    7.06e-8;    % C17-C18
    49.52,      1.630,   -6.35e-4,    7.98e-8;    % C19-C21
    57.03,      1.878,   -7.31e-4,    9.19e-8;    % C22-C24
    66.63,      2.194,   -8.55e-4,    1.07e-7;    % C25-C29
    82.18,      2.706,   -1.05e-3,    1.32e-7;    % C30-C37
   115.26,      3.795,   -1.48e-3,    1.86e-7     % C38-C80
]; 
H_ig_ref_per_mass = [-20; 20; 0; 7.5; 15; 17; 17; 25; 25; 33; 2; 2; 2; 2; 93; 93; 93; 31; 31; 31; 8; 8];
H_ig_ref = H_ig_ref_per_mass  .* M_gmol * R 

% H_ig_ref = [8330.8; 19459.1; 2.6; 9761.1; 19519.6; 29278.1; 29278.1; 39036.6; 39036.6; 48795.1; ...
%             58553.6; 68312.1; 78069.9; 86301.3; 105419.0; 131284.9; 158312.6; 180345.2; 209390.4; ...
%             246519.5; 306655.3; 434614.2]; 

% H_ig_per_MR = (-1342 + 83.67.*M_gmol) ; 
% H_ig_ref = H_ig_per_MR * R  ;  % [J/mol]


dTdh = 0.025;
h = [0;175;204;228;327]; 
 [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_nonisothermal(h(5), 175, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref)
% [comp_h, press_h, pressbub_h, pressdew_h] = main(h(5), 175, comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP, M_gmol, 0)

%%
[H_abs, H_mix, H_ig_spec, H_res_spec, H_abs_spec] = calculate_absolute_enthalpy(temp_ref, press_ref, comp_ref, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);

R = 8.314462618;
H_ig_plot = H_ig_spec  / R;
H_res_plot = H_res_spec / R;
H_tot_plot = H_abs_spec  / R;

% Plot
figure;
idx = 3:length(M_gmol); % Skip N2, CO2
plot(M_gmol(idx), H_res_plot(idx), 'b--', 'LineWidth', 2);
hold on;
plot (M_gmol(idx), H_ig_plot(idx), 'r-', 'LineWidth', 2);
plot (M_gmol(idx), H_tot_plot(idx), 'k-', 'LineWidth', 2.5);

xlim([0 700]);
ylim([-50 950]); 
xlabel('Molecular weight');
ylabel('H/M (K)');
legend('H^{res}', 'H^{ig}', 'H^{tot}', 'Location', 'northwest');
grid on;
