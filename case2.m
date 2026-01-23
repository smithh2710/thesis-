
components = {'C1', 'C2', 'C3', 'nC5', 'nC7', 'nC10'};

comp = [80.97; 5.66; 3.06; 4.57; 3.30; 2.44];
comp = comp / sum(comp);  % Normalize to mole fractions (sums to 1)

Pc = [46.00; 48.84; 42.46; 33.74; 27.36; 21.08] * 1e5;  % [Pa]
Tc = [190.6; 305.4; 369.8; 469.6; 540.2; 617.6];        % [K]
acentric = [0.008; 0.098; 0.152; 0.251; 0.351; 0.490];
M_gmol = [16.04; 30.07; 44.10; 72.15; 100.21; 142.29];
Vc = [99.2; 148.3; 203.0; 304.0; 432.0; 603.0];  % [cm³/mol]

n = length(comp);
R = 8.3144598;


BIP = zeros(n);
BIP(1, 2) = 0.0256;  % C1-C2
BIP(1, 3) = 0.0270;  % C1-C3
BIP(1, 4) = 0.0298;  % C1-nC5
BIP(1, 5) = 0.0326;  % C1-nC7
BIP(1, 6) = 0.0368;  % C1-nC10
BIP(2:end, 1) = BIP(1, 2:end)';


press_ref = 275e5;  % 275 bar
temp = 275;         % 275 K
h_ref = 1000;       % Reference depth [m]

depth_range = [1200 ,900 ]
vt_params.Vc = Vc;
vt_params.components = components;
[GOC_depth, GOC_pressure, GOC_comp, convergence_info] = detect_GOC(comp, press_ref, temp, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, 1, vt_params)

% [comp_h, press_h, pressbub_h, pressdew_h] = main(1006, h_ref, comp, press_ref, temp, Pc, Tc, acentric, BIP, M_gmol, 1, vt_params)
%%
vt_methods = [0, 1, 2, 3, 4, 5];
vt_names = {'No VT', 'Peneloux', 'Magoulas-Tassios', 'Ungerer-Batut', 'Baled', 'Abudour'};
n_methods = length(vt_methods);

% Colors for each method
colors = [0.0 0.0 0.0;    % Black - No VT
          0.0 0.4 0.8;    % Blue - Peneloux
          0.8 0.0 0.0;    % Red - MT
          0.0 0.6 0.3;    % Green - UB
          0.8 0.4 0.0;    % Orange - Baled
          0.6 0.0 0.6];   % Purple - Abudour

line_styles = {'-', '-', '--', '-.', '-', ':'};
line_widths = [2, 2, 1.5, 1.5, 2, 1.5];

%% DEPTH ARRAY
depths = 600:10:1200;
n_depths = length(depths);

% Convert to relative depth (h - h_ref)
% Positive = above reference (shallower)
% Negative = below reference (deeper)
depths_relative = h_ref - depths;  % h_ref - h, so deeper = negative

%% STORAGE
results = struct();
for m = 1:n_methods
    results(m).name = vt_names{m};
    results(m).P = zeros(n_depths, 1);
    results(m).comp = zeros(n_depths, n);
    results(m).rho = zeros(n_depths, 1);
end

%% CALCULATE FOR ALL METHODS
fprintf('=========================================================\n');
fprintf('CASE 2: Compositional Grading - All VT Methods\n');
fprintf('=========================================================\n');
fprintf('Reference: h = %d m, P = %.0f bar, T = %.0f K\n', h_ref, press_ref/1e5, temp);
fprintf('Depth range: %.0f - %.0f m\n', min(depths), max(depths));
fprintf('Relative depth: %.0f to %.0f m (ref = 0)\n', min(depths_relative), max(depths_relative));
fprintf('=========================================================\n\n');

for m = 1:n_methods
    vt = vt_methods(m);
    fprintf('Calculating %s (VT=%d)... ', vt_names{m}, vt);
    
    for i = 1:n_depths
        h = depths(i);
        
        try
            % Compositional grading
            [comp_h, P_h, ~, ~] = main(h, h_ref, comp, press_ref, temp, Pc, Tc, acentric, BIP, M_gmol, vt, vt_params);
            comp_h = comp_h(:);
            
            % Ensure composition is normalized (mole fractions sum to 1)
            comp_h = comp_h / sum(comp_h);
            
            results(m).P(i) = P_h;
            results(m).comp(i, :) = comp_h';
            
            % Calculate density
            results(m).rho(i) = calc_density(comp_h, P_h, temp, Pc, Tc, acentric, BIP, M_gmol, Vc, R, vt, vt_params);
        catch ME
            fprintf('\nError at depth %d: %s\n', h, ME.message);
            results(m).P(i) = NaN;
            results(m).comp(i, :) = NaN;
            results(m).rho(i) = NaN;
        end
    end
    
    fprintf('Done!\n');
end

%% VERIFY COMPOSITIONS ARE MOLE FRACTIONS
fprintf('\n--- Composition Verification ---\n');
for m = 1:n_methods
    comp_sum = sum(results(m).comp, 2);
    comp_max = max(results(m).comp(:));
    comp_min = min(results(m).comp(:));
    fprintf('%s: sum = %.4f-%.4f, range = [%.4f, %.4f]\n', ...
        vt_names{m}, min(comp_sum), max(comp_sum), comp_min, comp_max);
end

%% FIGURE 1: PRESSURE vs RELATIVE DEPTH (Book Style)
figure('Position', [50, 100, 600, 500], 'Color', 'w');
hold on;

for m = 1:n_methods
    plot(results(m).P/1e5, depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m), 'DisplayName', vt_names{m});
end

% Add reference line at h = 0
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

set(gca, 'FontSize', 11);
xlabel('Pressure [bar]', 'FontSize', 12);
ylabel('Height Relative to Reference [m]', 'FontSize', 12);
title('Case 2: Static Pressure vs Height', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
grid on;
box on;

% Add annotation
text(0.02, 0.98, sprintf('Reference: h_{ref} = %d m', h_ref), ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');

%% FIGURE 2: METHANE (C1) vs RELATIVE DEPTH (Book Style)
figure('Position', [100, 100, 600, 500], 'Color', 'w');
hold on;

for m = 1:n_methods
    plot(results(m).comp(:,1), depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m), 'DisplayName', vt_names{m});
end

% Add reference line at h = 0
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

set(gca, 'FontSize', 11);
xlabel('Methane (C_1) Mole Fraction [-]', 'FontSize', 12);
ylabel('Height Relative to Reference [m]', 'FontSize', 12);
title('Case 2: Methane Distribution vs Height', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 9);
xlim([0 1]);
grid on;
box on;

%% FIGURE 3: COMBINED  C2, nC10 vs RELATIVE DEPTH (Book Style)
figure('Position', [150, 100, 800, 600], 'Color', 'w');
hold on;

% Component indices
idx_C1 = 1;   % Methane
idx_C2 = 2;   % Ethane
idx_nC10 = 6; % n-Decane

% Plot each method with different colors, each component with different line style
for m = 1:n_methods
    % % C1 - Solid line
    % plot(results(m).comp(:, idx_C1), depths_relative, '-', 'Color', colors(m,:), ...
    %     'LineWidth', line_widths(m));
    
    % C2 - Dashed line
    plot(results(m).comp(:, idx_C2), depths_relative, '--', 'Color', colors(m,:), ...
        'LineWidth', line_widths(m));
    
    % nC10 - Dotted line
    plot(results(m).comp(:, idx_nC10), depths_relative, ':', 'Color', colors(m,:), ...
        'LineWidth', line_widths(m) + 0.5);
end

% Add reference line at h = 0
plot([0 1], [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

set(gca, 'FontSize', 11);
xlabel('Mole Fraction [-]', 'FontSize', 12);
ylabel('Height Relative to Reference [m]', 'FontSize', 12);
title('Case 2: Composition vs Height (C_1, C_2, nC_{10})', 'FontSize', 13);
xlim([0 1]);
grid on;
box on;

% Create custom legend
h_methods = zeros(n_methods, 1);
for m = 1:n_methods
    h_methods(m) = plot(NaN, NaN, '-', 'Color', colors(m,:), 'LineWidth', 2);
end

h_comp = zeros(3, 1);
h_comp(1) = plot(NaN, NaN, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
h_comp(2) = plot(NaN, NaN, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
h_comp(3) = plot(NaN, NaN, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 2.5);

legend([h_methods; h_comp], [vt_names, {'C_1', 'C_2', 'nC_{10}'}], ...
    'Location', 'eastoutside', 'FontSize', 9);

%% FIGURE 4: ETHANE (C2) vs RELATIVE DEPTH (Book Style)
figure('Position', [200, 100, 600, 500], 'Color', 'w');
hold on;

for m = 1:n_methods
    plot(results(m).comp(:,2), depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m), 'DisplayName', vt_names{m});
end

% Add reference line at h = 0
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

set(gca, 'FontSize', 11);
xlabel('Ethane (C_2) Mole Fraction [-]', 'FontSize', 12);
ylabel('Height Relative to Reference [m]', 'FontSize', 12);
title('Case 2: Ethane Distribution vs Height', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 9);
xlim([0 0.15]);
grid on;
box on;

%% FIGURE 5: n-DECANE (nC10) vs RELATIVE DEPTH (Book Style)
figure('Position', [250, 100, 600, 500], 'Color', 'w');
hold on;

for m = 1:n_methods
    plot(results(m).comp(:,6), depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m), 'DisplayName', vt_names{m});
end

% Add reference line at h = 0
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

set(gca, 'FontSize', 11);
xlabel('n-Decane (nC_{10}) Mole Fraction [-]', 'FontSize', 12);
ylabel('Height Relative to Reference [m]', 'FontSize', 12);
title('Case 2: n-Decane Distribution vs Height', 'FontSize', 13);
legend('Location', 'best', 'FontSize', 9);
xlim([0 0.15]);
grid on;
box on;

%% FIGURE 6: DENSITY vs RELATIVE DEPTH (Book Style)
figure('Position', [300, 100, 600, 500], 'Color', 'w');
hold on;

for m = 1:n_methods
    plot(results(m).rho, depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m), 'DisplayName', vt_names{m});
end

% Add reference line at h = 0
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

set(gca, 'FontSize', 11);
xlabel('Density [kg/m³]', 'FontSize', 12);
ylabel('Height Relative to Reference [m]', 'FontSize', 12);
title('Case 2: Density vs Height', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
grid on;
box on;

%% FIGURE 7: COMBINED 2x2 SUBPLOT (Book Style)
figure('Position', [350, 50, 1000, 800], 'Color', 'w');

% Pressure
subplot(2,2,1);
hold on;
for m = 1:n_methods
    plot(results(m).P/1e5, depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m));
end
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5);
set(gca, 'FontSize', 10);
xlabel('Pressure [bar]', 'FontSize', 11);
ylabel('Height Relative to Reference [m]', 'FontSize', 11);
title('(a) Pressure', 'FontSize', 12);
grid on;

% Composition (C1, C2, nC10 combined)
subplot(2,2,2);
hold on;
for m = 1:n_methods
    plot(results(m).comp(:,1), depths_relative, '-', 'Color', colors(m,:), 'LineWidth', line_widths(m));
    plot(results(m).comp(:,2), depths_relative, '--', 'Color', colors(m,:), 'LineWidth', line_widths(m));
    plot(results(m).comp(:,6), depths_relative, ':', 'Color', colors(m,:), 'LineWidth', line_widths(m)+0.5);
end
plot([0 1], [0 0], 'k--', 'LineWidth', 0.5);
set(gca, 'FontSize', 10);
xlabel('Mole Fraction [-]', 'FontSize', 11);
ylabel('Height Relative to Reference [m]', 'FontSize', 11);
title('(b) Composition (—C_1, --C_2, :nC_{10})', 'FontSize', 12);
xlim([0 1]);
grid on;

% Ethane & n-Decane (zoomed)
subplot(2,2,3);
hold on;
for m = 1:n_methods
    plot(results(m).comp(:,2), depths_relative, '--', 'Color', colors(m,:), 'LineWidth', line_widths(m));
    plot(results(m).comp(:,6), depths_relative, ':', 'Color', colors(m,:), 'LineWidth', line_widths(m)+0.5);
end
plot([0 0.12], [0 0], 'k--', 'LineWidth', 0.5);
set(gca, 'FontSize', 10);
xlabel('Mole Fraction [-]', 'FontSize', 11);
ylabel('Height Relative to Reference [m]', 'FontSize', 11);
title('(c) C_2 (--) and nC_{10} (:) Zoomed', 'FontSize', 12);
xlim([0 0.12]);
grid on;

% Density
subplot(2,2,4);
hold on;
for m = 1:n_methods
    plot(results(m).rho, depths_relative, line_styles{m}, 'Color', colors(m,:), ...
        'LineWidth', line_widths(m));
end
plot(xlim, [0 0], 'k--', 'LineWidth', 0.5);
set(gca, 'FontSize', 10);
xlabel('Density [kg/m³]', 'FontSize', 11);
ylabel('Height Relative to Reference [m]', 'FontSize', 11);
title('(d) Density', 'FontSize', 12);
grid on;

% Add legend to figure
h_leg = zeros(n_methods, 1);
for m = 1:n_methods
    h_leg(m) = plot(NaN, NaN, '-', 'Color', colors(m,:), 'LineWidth', line_widths(m));
end
legend(h_leg, vt_names, 'Position', [0.42, 0.01, 0.18, 0.05], 'FontSize', 8, 'NumColumns', 3);

sgtitle('Case 2: Compositional Grading - Effect of Volume Translation', 'FontSize', 14, 'FontWeight', 'bold');

%% PRINT SUMMARY TABLE - MOLE FRACTIONS
fprintf('\n=========================================================\n');
fprintf('SUMMARY: Values at Reference (h - h_{ref} = 0 m)\n');
fprintf('=========================================================\n');
fprintf('%-20s %10s %10s %10s %10s %12s\n', 'VT Method', 'P [bar]', 'x_C1', 'x_C2', 'x_nC10', 'rho [kg/m³]');
fprintf('%s\n', repmat('-', 1, 75));

[~, ref_idx] = min(abs(depths - h_ref));
for m = 1:n_methods
    fprintf('%-20s %10.2f %10.4f %10.4f %10.4f %12.2f\n', ...
        vt_names{m}, ...
        results(m).P(ref_idx)/1e5, ...
        results(m).comp(ref_idx,1), ...
        results(m).comp(ref_idx,2), ...
        results(m).comp(ref_idx,6), ...
        results(m).rho(ref_idx));
end

fprintf('\n=========================================================\n');
fprintf('SUMMARY: Values at Bottom (h - h_{ref} = %.0f m)\n', min(depths_relative));
fprintf('=========================================================\n');
fprintf('%-20s %10s %10s %10s %10s %12s\n', 'VT Method', 'P [bar]', 'x_C1', 'x_C2', 'x_nC10', 'rho [kg/m³]');
fprintf('%s\n', repmat('-', 1, 75));

for m = 1:n_methods
    fprintf('%-20s %10.2f %10.4f %10.4f %10.4f %12.2f\n', ...
        vt_names{m}, ...
        results(m).P(end)/1e5, ...
        results(m).comp(end,1), ...
        results(m).comp(end,2), ...
        results(m).comp(end,6), ...
        results(m).rho(end));
end

fprintf('\n=========================================================\n');
fprintf('SUMMARY: Values at Top (h - h_{ref} = %.0f m)\n', max(depths_relative));
fprintf('=========================================================\n');
fprintf('%-20s %10s %10s %10s %10s %12s\n', 'VT Method', 'P [bar]', 'x_C1', 'x_C2', 'x_nC10', 'rho [kg/m³]');
fprintf('%s\n', repmat('-', 1, 75));

for m = 1:n_methods
    fprintf('%-20s %10.2f %10.4f %10.4f %10.4f %12.2f\n', ...
        vt_names{m}, ...
        results(m).P(1)/1e5, ...
        results(m).comp(1,1), ...
        results(m).comp(1,2), ...
        results(m).comp(1,6), ...
        results(m).rho(1));
end

fprintf('=========================================================\n');


%% ========== HELPER FUNCTION ==========

function rho = calc_density(comp, P, T, Pc, Tc, acentric, BIP, M_gmol, Vc, R, vt_method, vt_params)
% Calculate density with volume translation
% CORRECTED: Abudour sign convention (no double negation)

    comp = comp(:);
    
    % Get Z-factor
    [~, Z] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP);
    if length(Z) > 1
        Z = min(Z);  % Liquid-like root
    end
    
    % Molar volume from EOS
    V_eos = Z * R * T / P;
    
    % Get c_mix from VT method
    c_mix = 0;
    switch vt_method
        case 0  % No VT
            c_mix = 0;
        case 1  % Peneloux
            [~, ~, mix] = peneloux_volume_shift(Pc, Tc, acentric, comp, Vc);
            c_mix = mix.c_mix;
        case 2  % Magoulas-Tassios
            [~, ~, mix] = magoulas_tassios_volume_shift(T, Pc, Tc, acentric, comp, Vc);
            c_mix = mix.c_mix;
        case 3  % Ungerer-Batut
            [~, ~, mix] = ungerer_batut_volume_shift(T, Pc, Tc, acentric, M_gmol, comp, Vc);
            c_mix = mix.c_mix;
        case 4  % Baled
            [~, ~, mix] = baled_volume_shift(T, Pc, Tc, acentric, M_gmol, vt_params.components, comp, Vc);
            c_mix = mix.c_mix;
        case 5  % Abudour - CORRECTED: No negation needed
            [~, ~, mix] = abudour_volume_shift(comp, P, T, Pc, Tc, acentric, Vc, vt_params.components, BIP, false);
            c_mix = mix.c_mix;  % Already in Peneloux convention
    end
    
    % Corrected volume
    V_corrected = V_eos - c_mix;
    
    % Density
    MW_mix = sum(comp .* M_gmol) / 1000;  % kg/mol
    rho = MW_mix / V_corrected;
end