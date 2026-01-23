function [GOC_depth, GOC_type, GOC_pressure, GOC_comp, info] = detect_UnGOC(comp_ref, press_ref, temp_ref, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params)
% DETECT_GOC_NEW - Detect Gas-Oil Contact depth and type
% =========================================================================
% This function implements the GOC detection algorithm based on:
%   - Hoier & Whitson (2000) methodology
%   - Saturation type tracking (Bubble vs Dew point dominated)
%   - Stability analysis for phase split detection
%
% GOC TYPES:
%   'saturated'     : P crosses Psat (phase split occurs)
%   'undersaturated': Saturation type changes (Pb = Pd, critical transition)
%   'none'          : No GOC found in depth range
%
% =========================================================================
% INPUTS:
%   comp_ref    : Reference composition [mole fractions] (column vector)
%   press_ref   : Reference pressure [Pa]
%   temp_ref    : Reference temperature [K]
%   h_ref       : Reference depth [m]
%   depth_range : [h_start, h_end] - Depth range to scan [m]
%   Pc          : Critical pressures [Pa]
%   Tc          : Critical temperatures [K]
%   acentric    : Acentric factors [-]
%   BIP         : Binary interaction parameter matrix
%   M_gmol      : Molecular weights [g/mol]
%   vt_type     : Volume translation method (0-5)
%   vt_params   : VT parameters structure
%
% OUTPUTS:
%   GOC_depth   : Depth of GOC [m] (NaN if not found)
%   GOC_type    : 'saturated', 'undersaturated', or 'none'
%   GOC_pressure: Pressure at GOC [Pa]
%   GOC_comp    : Composition at GOC [mole fractions]
%   info        : Structure with diagnostic information
%
% =========================================================================
% USAGE:
%   [GOC_depth, GOC_type, GOC_P, GOC_z, info] = detect_GOC_new(...
%       comp, P_ref, T_ref, h_ref, [h_start, h_end], ...
%       Pc, Tc, omega, BIP, M, vt_type, vt_params);
%
% =========================================================================

%% Initialize
comp_ref = comp_ref(:);
n = length(comp_ref);

% Default outputs
GOC_depth = NaN;
GOC_type = 'none';
GOC_pressure = NaN;
GOC_comp = comp_ref;

% Tolerances
tol_sat = 1e-8;         % Tolerance for saturation pressure calculation
maxiter_sat = 200;      % Max iterations for Psat
tol_GOC = 0.5;          % Depth tolerance for GOC refinement [m]

%% Parse depth range
h_start = min(depth_range);
h_end = max(depth_range);

% Determine scan direction (shallow to deep or deep to shallow)
if h_ref <= (h_start + h_end) / 2
    % Reference is near top, scan downward
    dh = 5;  % Step size [m]
    depths = h_ref:dh:h_end;
else
    % Reference is near bottom, scan upward
    dh = -5;
    depths = h_ref:dh:h_start;
end

n_depths = length(depths);

%% Step 1: Determine initial saturation type at reference
fprintf('=========================================================\n');
fprintf('GOC DETECTION ALGORITHM\n');
fprintf('=========================================================\n');
fprintf('Reference: h = %.1f m, P = %.2f bar, T = %.1f K\n', h_ref, press_ref/1e5, temp_ref);
fprintf('Scan range: %.1f to %.1f m\n', depths(1), depths(end));
fprintf('=========================================================\n\n');

[Pb_ref, Pd_ref, sat_type_ref] = get_saturation_type(comp_ref, temp_ref, Pc, Tc, acentric, BIP, press_ref, tol_sat, maxiter_sat);

fprintf('Initial saturation type: %s\n', sat_type_ref);
fprintf('  Pb = %.2f bar, Pd = %.2f bar, P = %.2f bar\n\n', Pb_ref/1e5, Pd_ref/1e5, press_ref/1e5);

%% Step 2: Scan depths
fprintf('Scanning depths...\n');
fprintf('%-8s %-10s %-10s %-10s %-10s %-10s %-12s\n', ...
    'Depth', 'P [bar]', 'Pb [bar]', 'Pd [bar]', 'Psat', 'Stable?', 'Sat Type');
fprintf('%s\n', repmat('-', 1, 75));

% Storage
info.depths = depths;
info.P = zeros(n_depths, 1);
info.Pb = zeros(n_depths, 1);
info.Pd = zeros(n_depths, 1);
info.sat_type = cell(n_depths, 1);
info.stable = true(n_depths, 1);
info.comp = zeros(n_depths, n);

prev_sat_type = sat_type_ref;
GOC_found = false;
GOC_bracket = [];  % [h_before, h_after] for refinement

for i = 1:n_depths
    h = depths(i);
    
    % Skip reference depth (already computed)
    if abs(h - h_ref) < 0.1
        info.P(i) = press_ref;
        info.Pb(i) = Pb_ref;
        info.Pd(i) = Pd_ref;
        info.sat_type{i} = sat_type_ref;
        info.stable(i) = true;
        info.comp(i, :) = comp_ref';
        
        fprintf('%-8.1f %-10.2f %-10.2f %-10.2f %-10.2f %-10s %-12s (ref)\n', ...
            h, press_ref/1e5, Pb_ref/1e5, Pd_ref/1e5, max(Pb_ref, Pd_ref)/1e5, 'Yes', sat_type_ref);
        continue;
    end
    
    % Get composition and pressure at depth h using CG calculation
    try
        [comp_h, P_h, Pb_h, Pd_h] = main(h, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
        comp_h = comp_h(:);
    catch ME
        fprintf('Warning: CG calculation failed at h = %.1f m: %s\n', h, ME.message);
        continue;
    end
    
    % Store results
    info.P(i) = P_h;
    info.comp(i, :) = comp_h';
    
    % Calculate saturation pressures if not returned properly
    if isnan(Pb_h) || Pb_h <= 0
        [Pb_h, ~] = calc_bubble_pressure(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    end
    if isnan(Pd_h) || Pd_h <= 0
        [Pd_h, ~] = calc_dew_pressure(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    end
    
    info.Pb(i) = Pb_h;
    info.Pd(i) = Pd_h;
    
    % Determine saturation type
    [~, ~, sat_type_h] = get_saturation_type(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    info.sat_type{i} = sat_type_h;
    
    % Determine Psat (the relevant saturation pressure)
    if strcmp(sat_type_h, 'Bubble')
        Psat = Pb_h;
    else
        Psat = Pd_h;
    end
    
    % Check stability (simple check: is P > Psat?)
    % More rigorous: use stability_analysis_ssi
    is_stable = (P_h > Psat) || isnan(Psat);
    info.stable(i) = is_stable;
    
    % Print status
    stable_str = 'Yes';
    if ~is_stable
        stable_str = 'NO';
    end
    fprintf('%-8.1f %-10.2f %-10.2f %-10.2f %-10.2f %-10s %-12s', ...
        h, P_h/1e5, Pb_h/1e5, Pd_h/1e5, Psat/1e5, stable_str, sat_type_h);
    
    %% Check for GOC conditions
    
    % Condition 1: SATURATED GOC - Phase becomes unstable (P crosses Psat)
    if ~is_stable && ~GOC_found
        fprintf(' *** UNSTABLE - Saturated GOC candidate\n');
        GOC_found = true;
        GOC_type = 'saturated';
        
        % Bracket for refinement
        if i > 1
            GOC_bracket = [depths(i-1), h];
        else
            GOC_bracket = [h - abs(dh), h];
        end
        break;  % Exit scan, proceed to refinement
    end
    
    % Condition 2: UNDERSATURATED GOC - Saturation type changed while still stable
    if is_stable && ~strcmp(sat_type_h, prev_sat_type) && ~GOC_found
        fprintf(' *** TYPE CHANGED - Undersaturated GOC candidate\n');
        GOC_found = true;
        GOC_type = 'undersaturated';
        
        % Bracket for refinement
        if i > 1
            GOC_bracket = [depths(i-1), h];
        else
            GOC_bracket = [h - abs(dh), h];
        end
        break;  % Exit scan, proceed to refinement
    end
    
    fprintf('\n');
    prev_sat_type = sat_type_h;
end

%% Step 3: Refine GOC location
if GOC_found && ~isempty(GOC_bracket)
    fprintf('\n=========================================================\n');
    fprintf('Refining GOC location...\n');
    fprintf('Bracket: [%.1f, %.1f] m\n', GOC_bracket(1), GOC_bracket(2));
    fprintf('=========================================================\n');
    
    if strcmp(GOC_type, 'saturated')
        % Use interval halving (bisection) to find where P = Psat
        [GOC_depth, GOC_pressure, GOC_comp] = refine_saturated_GOC_bisection(...
            GOC_bracket, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, ...
            tol_sat, maxiter_sat, tol_GOC);
        
    elseif strcmp(GOC_type, 'undersaturated')
        % Use secant method to find where Pb = Pd
        [GOC_depth, GOC_pressure, GOC_comp] = refine_undersaturated_GOC_secant(...
            GOC_bracket, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, ...
            tol_sat, maxiter_sat, tol_GOC);
    end
end

%% Output results
fprintf('\n=========================================================\n');
fprintf('GOC DETECTION RESULT\n');
fprintf('=========================================================\n');

if GOC_found
    fprintf('GOC Type:     %s\n', upper(GOC_type));
    fprintf('GOC Depth:    %.2f m\n', GOC_depth);
    fprintf('GOC Pressure: %.2f bar\n', GOC_pressure/1e5);
else
    fprintf('No GOC detected in range [%.1f, %.1f] m\n', h_start, h_end);
end
fprintf('=========================================================\n');

% Store final info
info.GOC_depth = GOC_depth;
info.GOC_type = GOC_type;
info.GOC_pressure = GOC_pressure;
info.GOC_comp = GOC_comp;
info.GOC_bracket = GOC_bracket;

end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [Pb, Pd, sat_type] = get_saturation_type(comp, T, Pc, Tc, acentric, BIP, P_init, tol, maxiter)
% GET_SATURATION_TYPE - Determine if fluid is bubble point or dew point type
%
% Returns 'Bubble' if Pb > Pd, 'Dew' otherwise

    comp = comp(:);
    
    % Calculate bubble point
    [Pb, ~] = calc_bubble_pressure(comp, T, Pc, Tc, acentric, BIP, P_init, tol, maxiter);
    
    % Calculate dew point
    [Pd, ~] = calc_dew_pressure(comp, T, Pc, Tc, acentric, BIP, P_init, tol, maxiter);
    
    % Determine type
    if ~isnan(Pb) && ~isnan(Pd)
        if Pb > Pd
            sat_type = 'Bubble';
        else
            sat_type = 'Dew';
        end
    elseif ~isnan(Pb)
        sat_type = 'Bubble';
    elseif ~isnan(Pd)
        sat_type = 'Dew';
    else
        sat_type = 'Unknown';
    end
end


function [Pb, y] = calc_bubble_pressure(comp, T, Pc, Tc, acentric, BIP, P_init, tol, maxiter)
% CALC_BUBBLE_PRESSURE - Calculate bubble point pressure with error handling

    try
        Pb_est = pressbubest_multicomp(comp, T, Pc, Tc, acentric);
        if isnan(Pb_est) || Pb_est <= 0
            Pb_est = P_init;
        end
        [Pb, y] = pressbub_multicomp_newton(comp, Pb_est, T, Pc, Tc, acentric, BIP, tol, maxiter);
        
        % Validate result
        if ~isreal(Pb) || Pb <= 0 || ~isfinite(Pb) || Pb > 1000e5
            Pb = NaN;
            y = comp;
        end
    catch
        Pb = NaN;
        y = comp;
    end
end


function [Pd, x] = calc_dew_pressure(comp, T, Pc, Tc, acentric, BIP, P_init, tol, maxiter)
% CALC_DEW_PRESSURE - Calculate dew point pressure with error handling

    try
        Pd_est = pressdewest_multicomp(comp, T, Pc, Tc, acentric);
        if isnan(Pd_est) || Pd_est <= 0
            Pd_est = P_init;
        end
        [Pd, x] = pressdew_multicomp_newton(comp, Pd_est, T, Pc, Tc, acentric, BIP, tol, maxiter);
        
        % Validate result
        if ~isreal(Pd) || Pd <= 0 || ~isfinite(Pd) || Pd > 1000e5
            Pd = NaN;
            x = comp;
        end
    catch
        Pd = NaN;
        x = comp;
    end
end


function [GOC_h, GOC_P, GOC_comp] = refine_saturated_GOC_bisection(bracket, h_ref, comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat, tol_GOC)
% REFINE_SATURATED_GOC_BISECTION - Find saturated GOC using interval halving
%
% Find depth where P = Psat (Delta_P = 0)
% Delta_P = (P - Psat) / P

    h_low = bracket(1);
    h_high = bracket(2);
    
    max_iter = 30;
    
    % Evaluate Delta_P at boundaries
    [Delta_P_low, ~, ~, ~] = calc_delta_P(h_low, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat);
    
    fprintf('Bisection: h_low = %.2f m, Delta_P = %.6f\n', h_low, Delta_P_low);
    
    GOC_h = (h_low + h_high) / 2;
    GOC_P = press_ref;
    GOC_comp = comp_ref;
    
    for iter = 1:max_iter
        h_mid = (h_low + h_high) / 2;
        
        [Delta_P_mid, P_mid, comp_mid, Psat_mid] = calc_delta_P(h_mid, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat);
        
        fprintf('  Iter %2d: h = %.2f m, P = %.2f bar, Psat = %.2f bar, Delta_P = %.6f\n', ...
            iter, h_mid, P_mid/1e5, Psat_mid/1e5, Delta_P_mid);
        
        % Check convergence
        if abs(Delta_P_mid) < 1e-4 || abs(h_high - h_low) < tol_GOC
            GOC_h = h_mid;
            GOC_P = P_mid;
            GOC_comp = comp_mid;
            fprintf('Converged: GOC at h = %.2f m\n', GOC_h);
            break;
        end
        
        % Bisection update
        if Delta_P_low * Delta_P_mid < 0
            h_high = h_mid;
        else
            h_low = h_mid;
            Delta_P_low = Delta_P_mid;
        end
        
        GOC_h = h_mid;
        GOC_P = P_mid;
        GOC_comp = comp_mid;
    end
end


function [GOC_h, GOC_P, GOC_comp] = refine_undersaturated_GOC_secant(bracket, h_ref, comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat, tol_GOC)
% REFINE_UNDERSATURATED_GOC_SECANT - Find undersaturated GOC using secant method
%
% Find depth where Pb = Pd (Delta_Psat = 0)
% Delta_Psat = Pb - Pd

    h0 = bracket(1);
    h1 = bracket(2);
    
    max_iter = 20;
    
    % Evaluate Delta_Psat at initial points
    [f0, ~, ~] = calc_delta_Psat(h0, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat);
    [f1, P1, comp1] = calc_delta_Psat(h1, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat);
    
    fprintf('Secant: h0 = %.2f m (Pb-Pd = %.2f bar), h1 = %.2f m (Pb-Pd = %.2f bar)\n', ...
        h0, f0/1e5, h1, f1/1e5);
    
    GOC_h = h1;
    GOC_P = P1;
    GOC_comp = comp1;
    
    for iter = 1:max_iter
        % Secant update
        if abs(f1 - f0) < 1e-10
            fprintf('Warning: Secant denominator too small\n');
            break;
        end
        
        h_new = h1 - f1 * (h1 - h0) / (f1 - f0);
        
        % Keep within bounds
        h_new = max(min(bracket), min(max(bracket), h_new));
        
        [f_new, P_new, comp_new] = calc_delta_Psat(h_new, h_ref, comp_ref, press_ref, temp_ref, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat);
        
        fprintf('  Iter %2d: h = %.2f m, Pb-Pd = %.4f bar\n', iter, h_new, f_new/1e5);
        
        % Check convergence
        if abs(f_new) < 0.01e5 || abs(h_new - h1) < tol_GOC  % 0.01 bar tolerance
            GOC_h = h_new;
            GOC_P = P_new;
            GOC_comp = comp_new;
            fprintf('Converged: Undersaturated GOC at h = %.2f m\n', GOC_h);
            break;
        end
        
        % Update for next iteration
        h0 = h1;
        f0 = f1;
        h1 = h_new;
        f1 = f_new;
        
        GOC_h = h_new;
        GOC_P = P_new;
        GOC_comp = comp_new;
    end
end


function [Delta_P, P_h, comp_h, Psat] = calc_delta_P(h, h_ref, comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat)
% CALC_DELTA_P - Calculate normalized pressure difference from saturation
%
% Delta_P = (P - Psat) / P
% Positive: undersaturated (single phase)
% Negative: should be two-phase (unstable)

    % Get composition and pressure at depth h
    [comp_h, P_h, Pb_h, Pd_h] = main(h, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    comp_h = comp_h(:);
    
    % Get saturation pressures if needed
    if isnan(Pb_h) || Pb_h <= 0
        [Pb_h, ~] = calc_bubble_pressure(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    end
    if isnan(Pd_h) || Pd_h <= 0
        [Pd_h, ~] = calc_dew_pressure(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    end
    
    % Use the higher saturation pressure
    if ~isnan(Pb_h) && ~isnan(Pd_h)
        Psat = max(Pb_h, Pd_h);
    elseif ~isnan(Pb_h)
        Psat = Pb_h;
    elseif ~isnan(Pd_h)
        Psat = Pd_h;
    else
        Psat = P_h;  % Fallback
    end
    
    % Calculate Delta_P
    Delta_P = (P_h - Psat) / P_h;
end


function [Delta_Psat, P_h, comp_h] = calc_delta_Psat(h, h_ref, comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, tol_sat, maxiter_sat)
% CALC_DELTA_PSAT - Calculate difference between bubble and dew pressure
%
% Delta_Psat = Pb - Pd
% Zero at critical point (undersaturated GOC)

    % Get composition and pressure at depth h
    [comp_h, P_h, ~, ~] = main(h, h_ref, comp_ref, press_ref, temp_ref, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    comp_h = comp_h(:);
    
    % Calculate both saturation pressures
    [Pb_h, ~] = calc_bubble_pressure(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    [Pd_h, ~] = calc_dew_pressure(comp_h, temp_ref, Pc, Tc, acentric, BIP, P_h, tol_sat, maxiter_sat);
    
    % Handle NaN cases
    if isnan(Pb_h)
        Pb_h = 0;
    end
    if isnan(Pd_h)
        Pd_h = 0;
    end
    
    Delta_Psat = Pb_h - Pd_h;
end