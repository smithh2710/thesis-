% components = {'C2' , 'C3'}
% comp = [0.3094 ; 0.6906 ] ;
% press = 14.20320002 * 10^(5) ; % pa 
% temp = 288.75 ; % k 
% R = 8.3144598 ; 
% Pc = [48.84 ; 42.46]* 10^(5);
% Tc = [305.322 ; 369.89];
% acentric = [0.0995 ; 0.1521]; 
% vc = [145.8387 ; 200 ]; 
% M = [ 30.069 ; 44.0956 ] ; 

% BIP = zeros(2);
% % [s, c, mix] = baled_volume_shift(temp, Pc, Tc, acentric, M, components, comp, vc)  
% 
% % [s, c, mix] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol, comp, vc)
% 
% % [s, c, mix] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric, comp, vc)
% 
% % [s, c, mix] = peneloux_volume_shift(Pc, Tc, acentric, comp, vc)   
 
% [~, Z] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP);

% v_pr = ( Z * R * temp / press ) * 10^6 % [cm³/mol]
% vt_pr = v_pr - (mix.c_mix_cm3)




%% Demo for pressbub_multicomp_newton() and pressdew_multicomp_newton()
% CORRECTED VERSION:
%   - Abudour sign convention fixed (no double negation)
%   - Mole fractions displayed as 0-1 (not percentages)

components = {'CO2','N2','C1','C2','C3','iC4','nC4','iC5','nC5','nC6','C11','C33'};
comp = [0.12; 0.13; 54.98; 5.44; 4.88; 1.88; 2.48; 1.45; 1.30; 1.91; 21.53; 3.90];
comp = comp ./ sum(comp);  % Normalize to mole fractions (0-1)

Pc = [73.76; 33.94; 46.00; 48.84; 42.46; 36.48; 38.00; 33.84; 33.74; 29.69; 23.43; 9.90]*1e5;
Tc = [304.2; 126.2; 190.6; 305.4; 369.8; 408.1; 425.2; 460.4; 469.6; 507.4; 642.51; 929.34];
acentric = [0.2250; 0.0400; 0.0080; 0.0980; 0.1520; 0.1760; 0.1930; 0.2270; 0.2510; 0.2960; 0.4882; 1.0854];
M_gmol = [44.01; 28.01; 16.04; 30.07; 44.10; 58.12; 58.12; 72.15; 72.15; 86.18; 151.16; 459.34];
Vc = [94; 89.80; 99; 148; 203; 263; 255; 306; 304; 370; 668.76; 1998.60];

n = length(comp);
R = 8.3144598;

BIP = zeros(n);
BIP(3,4:end) = [0.0256, 0.0270, 0.0284, 0.0284, 0.0298, 0.0298, 0.0312, 0.0486, 0.085];
BIP(4:end,3) = BIP(3,4:end);

% BIP(1,3:end)= [0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.1,0.1]; 
% BIP(3:end,1 )= BIP(1,3:end);
% 
% BIP(2,1) = -0.0170 ; BIP(1,2)= BIP(2,1);
% BIP(2,3:end)=[0.0311,0.0515,0.0852,0.1033,0.0800,0.0922,0.1000,0.0800,0.0800,0.0800];
% BIP(3:end,2)=BIP(2,3:end) ; 
vt_params.Vc = Vc;
vt_params.components = components

h_ref = 100;
press_ref = 293.02e5;
temp = 361;
depth_range = [40, 100];
 
[GOC_depth, GOC_pressure, GOC_comp, info] = detect_GOC(comp, press_ref, temp, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, 1, vt_params)

% [comp_h, press_h, pressbub_h, pressdew_h] = main(1006, h_ref, comp, press_ref, temp, Pc, Tc, acentric, BIP, M_gmol, 1, vt_params)

%  pressbub_ini_h = pressbubest_multicomp(comp, temp, Pc, Tc, acentric) ; 
% [pressbub_h, ~] = pressbub_multicomp_newton(comp, pressbub_ini_h, temp, Pc, Tc, acentric, BIP, 10^(-12) , 500)

%% VT METHODS
vt_methods = [0, 1, 2, 3, 4, 5];
vt_names = {'No VT', 'Peneloux', 'Magoulas-Tassios', 'Ungerer-Batut', 'Baled', 'Abudour'};
n_methods = length(vt_methods);

% Colors for each method
colors = [0.0 0.0 0.0;   % Black - No VT
          0.0 0.4 0.8;   % Blue - Peneloux
          0.9 0.5 0.0;   % Orange - Magoulas-Tassios
          0.5 0.0 0.5;   % Purple - Ungerer-Batut
          0.0 0.6 0.3;   % Green - Baled
          0.8 0.0 0.0];  % Red - Abudour

results = struct();

%% RUN FOR EACH VT METHOD
for m = 1:n_methods
    vt = vt_methods(m);
    fprintf('Running: %s...\n', vt_names{m});
    
    % Detect GOC
    evalc('[GOC_depth, ~, ~, ~, ~] = detect_GOC(comp, press_ref, temp, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, vt, vt_params);');
    
    % Height arrays
    h_liq = unique([h_ref:-1:ceil(GOC_depth), GOC_depth], 'stable')';
    h_gas = unique([GOC_depth, floor(GOC_depth):-1:depth_range(1)], 'stable')';
    
    n_liq = length(h_liq);
    n_gas = length(h_gas);
    
    % Preallocate
    P_liq = zeros(n_liq,1); P_bub = zeros(n_liq,1); x_liq = zeros(n_liq,n); rho_liq = zeros(n_liq,1);
    P_gas = zeros(n_gas,1); P_dew = zeros(n_gas,1); x_gas = zeros(n_gas,n); rho_gas = zeros(n_gas,1);
    
    % Liquid zone
    for i = 1:n_liq
        [x_liq(i,:), P_liq(i), P_bub(i), ~] = main(h_liq(i), h_ref, comp, press_ref, temp, Pc, Tc, acentric, BIP, M_gmol, vt, vt_params);
        
        % Normalize composition to ensure sum = 1
        x_liq(i,:) = x_liq(i,:) / sum(x_liq(i,:));
        
        rho_liq(i) = calc_density(x_liq(i,:)', P_liq(i), temp, Pc, Tc, acentric, BIP, M_gmol, Vc, R, vt, vt_params, 'liquid');
    end
    
    % Flash at GOC
    comp_goc = x_liq(end,:)';
    press_goc = P_liq(end);
    [~, comp_vap, ~, ~] = vaporliquideq(press_goc, temp, comp_goc, Pc, Tc, acentric, BIP, 1e-8, 150);
    
    % Gas zone
    for i = 1:n_gas
        [x_gas(i,:), P_gas(i), ~, P_dew(i)] = main(h_gas(i), GOC_depth, comp_vap, press_goc, temp, Pc, Tc, acentric, BIP, M_gmol, vt, vt_params);
        
        % Normalize composition to ensure sum = 1
        x_gas(i,:) = x_gas(i,:) / sum(x_gas(i,:));
        
        rho_gas(i) = calc_density(x_gas(i,:)', P_gas(i), temp, Pc, Tc, acentric, BIP, M_gmol, Vc, R, vt, vt_params, 'vapor');
    end
    
    % Combine
    results(m).heights = [h_liq; h_gas(2:end)];
    results(m).P = [P_liq; P_gas(2:end)];
    results(m).P_sat = [P_bub; P_dew(2:end)];
    results(m).x = [x_liq; x_gas(2:end,:)];
    results(m).rho = [rho_liq; rho_gas(2:end)];
    results(m).GOC = GOC_depth;
    results(m).name = vt_names{m};
end
fprintf('Done!\n\n');

%% VERIFY COMPOSITIONS ARE MOLE FRACTIONS (0-1)
fprintf('--- Composition Verification ---\n');
for m = 1:n_methods
    comp_sum = sum(results(m).x, 2);
    comp_max = max(results(m).x(:));
    comp_min = min(results(m).x(:));
    fprintf('%s: sum = %.4f-%.4f, range = [%.4f, %.4f]\n', ...
        vt_names{m}, min(comp_sum), max(comp_sum), comp_min, comp_max);
end
fprintf('\n');

%% PLOT 1: Pressure & Saturation Pressure vs Depth
figure('Position', [50, 100, 600, 500], 'Color', 'w'); hold on;

for m = 1:n_methods
    plot(results(m).P/1e5, results(m).heights, '-', 'Color', colors(m,:), 'LineWidth', 1.5);
    plot(results(m).P_sat/1e5, results(m).heights, '--', 'Color', colors(m,:), 'LineWidth', 1.2);
end

set(gca, 'YDir', 'reverse');
xlabel('Pressure [bar]'); ylabel('Depth [m]');
title('Pressure & Saturation Pressure vs Depth');

% Custom legend
h_leg = zeros(n_methods, 1);
for m = 1:n_methods
    h_leg(m) = plot(NaN, NaN, '-', 'Color', colors(m,:), 'LineWidth', 2);
end
legend(h_leg, vt_names, 'Location', 'southeast');
grid on;

% Add annotation for line styles
text(0.02, 0.02, 'Solid: P_{res}, Dashed: P_{sat}', 'Units', 'normalized', 'FontSize', 9);

%% PLOT 2: Composition (C1, C2, C7+) vs Depth - MOLE FRACTIONS (0-1)
figure('Position', [100, 100, 600, 500], 'Color', 'w'); hold on;

idx_C1 = 3; idx_C2 = 4; idx_C7plus = [11, 12];

for m = 1:n_methods
    x_C1 = results(m).x(:, idx_C1);           % Mole fraction (0-1)
    x_C2 = results(m).x(:, idx_C2);           % Mole fraction (0-1)
    x_C7plus = sum(results(m).x(:, idx_C7plus), 2);  % Mole fraction (0-1)
    
    plot(x_C1, results(m).heights, '-', 'Color', colors(m,:), 'LineWidth', 1.5);      % Solid: C1
    plot(x_C2, results(m).heights, '--', 'Color', colors(m,:), 'LineWidth', 1.2);     % Dashed: C2
    plot(x_C7plus, results(m).heights, ':', 'Color', colors(m,:), 'LineWidth', 1.5);  % Dotted: C7+
end

set(gca, 'YDir', 'reverse');
xlabel('Mole Fraction [-]'); ylabel('Depth [m]');
title('Composition vs Depth');
xlim([0 1]);  % Mole fraction range

% Custom legend for methods
h_leg = zeros(n_methods, 1);
for m = 1:n_methods
    h_leg(m) = plot(NaN, NaN, '-', 'Color', colors(m,:), 'LineWidth', 2);
end
legend(h_leg, vt_names, 'Location', 'best');
grid on;

% Add annotation for line styles
text(0.02, 0.02, 'Solid: C_1, Dashed: C_2, Dotted: C_{7+}', 'Units', 'normalized', 'FontSize', 9);

%% PLOT 3: Density vs Depth
figure('Position', [150, 100, 600, 500], 'Color', 'w'); hold on;

for m = 1:n_methods
    plot(results(m).rho, results(m).heights, '-', 'Color', colors(m,:), 'LineWidth', 1.5);
end

set(gca, 'YDir', 'reverse');
xlabel('Density [kg/m³]'); ylabel('Depth [m]');
title('Mixture Density vs Depth');
legend(vt_names, 'Location', 'southeast');
grid on;

%% SUMMARY TABLE
fprintf('=== GOC DEPTH COMPARISON ===\n');
fprintf('%-20s %10s\n', 'VT Method', 'GOC [m]');
fprintf('%s\n', repmat('-', 1, 32));
for m = 1:n_methods
    fprintf('%-20s %10.2f\n', vt_names{m}, results(m).GOC);
end

%% COMPOSITION SUMMARY (Mole Fractions)
fprintf('\n=== COMPOSITION AT REFERENCE DEPTH (h = %d m) ===\n', h_ref);
fprintf('%-20s %10s %10s %10s\n', 'VT Method', 'x_C1', 'x_C2', 'x_C7+');
fprintf('%s\n', repmat('-', 1, 52));
for m = 1:n_methods
    [~, idx] = min(abs(results(m).heights - h_ref));
    x_C1 = results(m).x(idx, 3);
    x_C2 = results(m).x(idx, 4);
    x_C7plus = sum(results(m).x(idx, [11,12]));
    fprintf('%-20s %10.4f %10.4f %10.4f\n', vt_names{m}, x_C1, x_C2, x_C7plus);
end


%% ========== HELPER FUNCTIONS ==========

function rho = calc_density(comp, P, T, Pc, Tc, acentric, BIP, M_gmol, Vc, R, vt_method, vt_params, phase)
% CALC_DENSITY - Calculate mixture density with proper VT mixing rules
%
% CORRECTED: Abudour sign convention (no double negation)
%
% Standard convention: V_corrected = V_EOS - c_mix (positive c reduces volume)

    % Ensure column vectors
    comp = comp(:);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    Vc = Vc(:);
    
    n = length(comp);
    
    % Get Z-factor for correct phase
    if strcmp(phase, 'liquid')
        [~, Z] = fugacitycoef_multicomp_liquid(comp, P, T, Pc, Tc, acentric, BIP);
    else
        [~, Z] = fugacitycoef_multicomp_vapor(comp, P, T, Pc, Tc, acentric, BIP);
    end
    
    % Ensure Z is scalar
    if length(Z) > 1
        if strcmp(phase, 'liquid')
            Z = min(Z);
        else
            Z = max(Z);
        end
    end
    
    % Molar volume from PR-EOS [m³/mol]
    V_eos = Z * R * T / P;
    
    % Get c_mix using proper mixing rule from VT function
    c_mix = get_cmix_from_vt(comp, P, T, Pc, Tc, acentric, M_gmol, Vc, vt_method, vt_params, BIP);
    
    % Corrected molar volume (SUBTRACT positive c for volume reduction)
    V_corrected = V_eos - c_mix;
    
    % Mixture molecular weight [kg/mol]
    MW_mix = sum(comp .* M_gmol) / 1000;
    
    % Density [kg/m³]
    rho = MW_mix / V_corrected;

end


function c_mix = get_cmix_from_vt(comp, P, T, Pc, Tc, acentric, M_gmol, Vc, vt_method, vt_params, BIP)
% GET_CMIX_FROM_VT - Get mixture volume translation using proper mixing rules

%   c_mix : Mixture volume translation [m³/mol] (positive = volume reduction)

    n = length(comp);
    
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
            components = vt_params.components;
            [~, ~, mix] = baled_volume_shift(T, Pc, Tc, acentric, M_gmol, components, comp, Vc);
            c_mix = mix.c_mix;
            
        case 5  % Abudour 
            components = vt_params.components;
            [~, ~, mix] = abudour_volume_shift(comp, P, T, Pc, Tc, acentric, Vc, components, BIP, false);
            c_mix = mix.c_mix;  
            
        otherwise
            c_mix = 0;
    end

end