%% RESERVOIR 1 - COMPLETE Isothermal Compositional Grading with GOR
% Reference: Pedersen et al. (2015) SPE-175085-MS
% =========================================================================
% This script generates COMPLETE profiles including:
%   - Oil zone (below GOC): Bubble point tracking
%   - Gas zone (above GOC): Dew point tracking via flash at GOC
%   - GOR calculation via flash to STO conditions
%   - Matches Figure 1 and Figure 3 format from the paper
% =========================================================================

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
M_gmol = [28.0; 44.0; 16.0; 30.1; 44.1; 58.1; 58.1; 72.2; 72.2; 86.2; ...
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
R = 8.3144598;
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
temp = 93 + 273.15;


vt_params.Vc = Vc;
vt_params.components = components;
vt_params.c_custom = c_JY;

fprintf('=========================================================================\n');
fprintf('STEP 1: GOC DETECTION\n');
fprintf('=========================================================================\n');

depth_range_goc = [0, 175];
[GOC_depth, GOC_pressure, GOC_comp, info] = detect_GOC(comp_ref, press_ref, temp,h_ref, depth_range_goc, Pc, Tc, acentric, BIP, M_gmol, 6, vt_params);

fprintf('\n*** GOC RESULT ***\n');
fprintf('Calculated GOC: %.2f m\n', GOC_depth);

%% ========================================================================
% STEP 2: FLASH AT GOC TO GET VAPOR COMPOSITION
% =========================================================================
fprintf('\n=========================================================================\n');
fprintf('STEP 2: FLASH CALCULATION AT GOC\n');
fprintf('=========================================================================\n');

[comp_at_goc, P_at_goc, ~, ~] = main(GOC_depth, h_ref, comp_ref, press_ref, temp, ...
    Pc, Tc, acentric, BIP, M_gmol, 6, vt_params);

try
    [K_goc, comp_vap_goc, comp_liq_goc, vapor_frac] = vaporliquideq(P_at_goc, temp, ...
        comp_at_goc, Pc, Tc, acentric, BIP, 1e-10, 200);

    fprintf('Flash at GOC successful.\n');
    fprintf('  Vapor C1: %.2f%%, Liquid C1: %.2f%%\n', comp_vap_goc(3)*100, comp_liq_goc(3)*100);

    comp_gas_ref = comp_vap_goc;
    press_gas_ref = P_at_goc;
    h_gas_ref = GOC_depth;
catch ME
    fprintf('Flash failed, using Wilson K-values.\n');
    K_wilson = wilsoneq(P_at_goc, temp, Pc, Tc, acentric);
    comp_gas_ref = K_wilson .* comp_at_goc;
    comp_gas_ref = comp_gas_ref / sum(comp_gas_ref);
    press_gas_ref = P_at_goc;
    h_gas_ref = GOC_depth;
end


fprintf('\n=========================================================================\n');
fprintf('STEP 3: OIL ZONE PROFILE\n');
fprintf('=========================================================================\n');

depths_oil = unique([GOC_depth, ceil(GOC_depth):5:350]);
n_oil = length(depths_oil);

oil_results.depth = zeros(n_oil, 1);
oil_results.P = zeros(n_oil, 1);
oil_results.Pbub = zeros(n_oil, 1);
oil_results.comp = zeros(n_oil, n);
oil_results.GOR = zeros(n_oil, 1);

fprintf('Calculating oil zone with GOR...\n');
for i = 1:n_oil
    h = depths_oil(i);

    [comp_h, P_h, Pb_h, ~] = main(h, h_ref, comp_ref, press_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, 6, vt_params);

    oil_results.depth(i) = h;
    oil_results.P(i) = P_h / 1e5;
    oil_results.Pbub(i) = Pb_h / 1e5;
    oil_results.comp(i, :) = comp_h(:)';

    % Calculate GOR via flash to STO
    try
        [GOR_h, ~, ~] = calculate_GOR_STO(comp_h, P_h, temp, Pc, Tc, acentric, BIP, M_gmol);
        oil_results.GOR(i) = GOR_h;
    catch
        oil_results.GOR(i) = NaN;
    end

    if mod(i, 10) == 0
        fprintf('  Depth %.1f m: P = %.1f bar, GOR = %.1f\n', h, P_h/1e5, oil_results.GOR(i));
    end
end
fprintf('Oil zone complete.\n');


fprintf('\n=========================================================================\n');
fprintf('STEP 4: GAS ZONE PROFILE\n');
fprintf('=========================================================================\n');

depths_gas = unique([0:5:floor(GOC_depth), GOC_depth]);
n_gas = length(depths_gas);

gas_results.depth = zeros(n_gas, 1);
gas_results.P = zeros(n_gas, 1);
gas_results.Pdew = zeros(n_gas, 1);
gas_results.comp = zeros(n_gas, n);
gas_results.GOR = zeros(n_gas, 1);

fprintf('Calculating gas zone with GOR...\n');
for i = 1:n_gas
    h = depths_gas(i);

    [comp_h, P_h, ~, Pd_h] = main(h, h_gas_ref, comp_gas_ref, press_gas_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, 6, vt_params);

    gas_results.depth(i) = h;
    gas_results.P(i) = P_h / 1e5;
    gas_results.Pdew(i) = Pd_h / 1e5;
    gas_results.comp(i, :) = comp_h(:)';

    % Calculate GOR via flash to STO
    try
        [GOR_h, ~, ~] = calculate_GOR_STO(comp_h, P_h, temp, Pc, Tc, acentric, BIP, M_gmol);
        gas_results.GOR(i) = GOR_h;
    catch
        gas_results.GOR(i) = NaN;
    end
end
fprintf('Gas zone complete.\n');



exp.depths = [0; 175; 204; 228; 327];
exp.P_bar = [279; 284; 286; 287; 293];
exp.Psat_bar = [270; 272; 267; 265; 242];
exp.GOR = [3739; 281.5; 280.9; 264.5; 227.8];

exp.C1 = [75.66; 50.04; 49.88; 48.89; 45.66];
exp.C7 = [1.29; 3.50; 3.49; 3.62; 4.01];
exp.C8 = [1.06; 3.75; 3.73; 3.85; 4.28];
exp.C9 = [0.54; 2.28; 2.31; 2.36; 2.58];
exp.C10plus = [1.57; 15.88; 16.11; 16.70; 17.66];
exp.C7plus = exp.C7 + exp.C8 + exp.C9 + exp.C10plus;

exp.N2 = [1.06; 0.42; 0.42; 0.39; 0.43];
exp.CO2 = [0.70; 0.69; 0.71; 0.70; 0.77];
exp.inorganics = exp.N2 + exp.CO2;

fprintf('\n=========================================================================\n');
fprintf('GENERATING PLOTS\n');
fprintf('=========================================================================\n');

% Combine data
all_depths = [gas_results.depth; oil_results.depth];
all_P = [gas_results.P; oil_results.P];
all_Psat = [gas_results.Pdew; oil_results.Pbub];
all_comp = [gas_results.comp; oil_results.comp];
all_GOR = [gas_results.GOR; oil_results.GOR];


figure('Position', [50, 100, 700, 550], 'Color', 'w');
hold on;

plot(all_P, all_depths, 'b-', 'LineWidth', 2, 'DisplayName', 'Sim Reservoir Pressure');
plot(all_Psat, all_depths, 'b--', 'LineWidth', 2, 'DisplayName', 'Sim Saturation Pressure');
plot(exp.P_bar, exp.depths, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'Exp Reservoir Pressure');
plot(exp.Psat_bar, exp.depths, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'Exp Saturation Pressure');

yline(GOC_depth, 'r:', 'LineWidth', 2, 'DisplayName', sprintf('Calc GOC: %.1f m', GOC_depth));
yline(140, 'r--', 'LineWidth', 2, 'DisplayName', 'Exp GOC: 140 m');

set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('Pressure [bar]', 'FontSize', 12);
ylabel('Depth [m]', 'FontSize', 12);
title('Reservoir #1: No T-gradient (Isothermal)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 9);
grid on;
xlim([220, 300]);
ylim([0, 350]);

figure('Position', [100, 100, 700, 550], 'Color', 'w');
hold on;

% Plot calculated GOR (log scale)
valid_GOR = all_GOR > 0 & isfinite(all_GOR);
semilogy(all_GOR(valid_GOR), all_depths(valid_GOR), 'b-', 'LineWidth', 2, ...
         'DisplayName', 'Sim GOR');

% Plot experimental GOR
semilogy(exp.GOR, exp.depths, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
         'DisplayName', 'Exp GOR');

% GOC lines
yline(GOC_depth, 'k:', 'LineWidth', 1.5, 'DisplayName', sprintf('Calc GOC: %.1f m', GOC_depth));
yline(140, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Exp GOC: 140 m');

set(gca, 'YDir', 'reverse', 'FontSize', 11);
xlabel('GOR [Sm³/Sm³]', 'FontSize', 12);
ylabel('Depth [m]', 'FontSize', 12);
title('Reservoir #1: GOR vs Depth (Isothermal)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 10);
grid on;
ylim([0, 350]);

figure('Position', [150, 50, 1200, 900], 'Color', 'w');

% Subplot 1: Pressure
subplot(2, 2, 1);
hold on;
plot(all_P, all_depths, 'b-', 'LineWidth', 2, 'DisplayName', 'P_{res} (calc)');
plot(all_Psat, all_depths, 'r--', 'LineWidth', 2, 'DisplayName', 'P_{sat} (calc)');
plot(exp.P_bar, exp.depths, 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b', ...
     'DisplayName', 'P_{res} (exp)');
plot(exp.Psat_bar, exp.depths, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'P_{sat} (exp)');
yline(GOC_depth, 'k:', 'LineWidth', 1.5);
yline(140, 'k--', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('Pressure [bar]', 'FontSize', 11);
ylabel('Depth [m]', 'FontSize', 11);
title('(a) Pressure vs Depth', 'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 8);
grid on;
xlim([220, 300]);

% Subplot 2: C1
subplot(2, 2, 2);
hold on;
plot(all_comp(:,3)*100, all_depths, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
plot(exp.C1, exp.depths, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'Experimental');
yline(GOC_depth, 'k:', 'LineWidth', 1.5);
yline(140, 'k--', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('C_1 Mole %', 'FontSize', 11);
ylabel('Depth [m]', 'FontSize', 11);
title('(b) Methane vs Depth', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9);
grid on;

% Subplot 3: C7+
subplot(2, 2, 3);
hold on;
C7plus_calc = sum(all_comp(:, 11:end), 2) * 100;
plot(C7plus_calc, all_depths, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
plot(exp.C7plus, exp.depths, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'Experimental');
yline(GOC_depth, 'k:', 'LineWidth', 1.5);
yline(140, 'k--', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('C_{7+} Mole %', 'FontSize', 11);
ylabel('Depth [m]', 'FontSize', 11);
title('(c) C_{7+} vs Depth', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9);
grid on;

% Subplot 4: GOR (log scale)
subplot(2, 2, 4);
hold on;
valid = all_GOR > 0 & isfinite(all_GOR);
semilogy(all_GOR(valid), all_depths(valid), 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated');
semilogy(exp.GOR, exp.depths, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
         'DisplayName', 'Experimental');
yline(GOC_depth, 'k:', 'LineWidth', 1.5);
yline(140, 'k--', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('GOR [Sm³/Sm³]', 'FontSize', 11);
ylabel('Depth [m]', 'FontSize', 11);
title('(d) GOR vs Depth', 'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 9);
grid on;

sgtitle('Reservoir 1 - Isothermal Compositional Grading (Pedersen 2015)', ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n=========================================================================\n');
fprintf('GOR COMPARISON\n');
fprintf('=========================================================================\n');
fprintf('%-10s %15s %15s %15s\n', 'Depth [m]', 'GOR_calc', 'GOR_exp', 'Ratio');
fprintf('%s\n', repmat('-', 1, 58));

for i = 1:length(exp.depths)
    d = exp.depths(i);

    % Find closest calculated depth
    [~, idx] = min(abs(all_depths - d));
    GOR_calc = all_GOR(idx);
    GOR_exp = exp.GOR(i);

    if GOR_exp > 0
        ratio = GOR_calc / GOR_exp;
    else
        ratio = NaN;
    end

    fprintf('%-10.0f %15.1f %15.1f %15.3f\n', d, GOR_calc, GOR_exp, ratio);
end

fprintf('\n=========================================================================\n');
fprintf('COMPLETE\n');
fprintf('=========================================================================\n');


