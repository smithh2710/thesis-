function [GOC_depth, GOC_pressure, GOC_comp, info] = detect_GOC(comp_ref, press_ref, temp, h_ref, depth_range, Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params)
% DETECT_GOC - Gas-Oil Contact Detection Algorithm (great for saturated GOC)
% =========================================================================
% Complete implementation based on Nikpoor (2014) Section 4.4 and Figure 4.2
% Handles BOTH saturated and undersaturated GOC cases
%
% GOC Types:
% ----------
% 1. SATURATED GOC: Reservoir pressure equals saturation pressure
%    - P_res = P_sat (either P_bub or P_dew)
%    - Detection: Sign change in ΔP = (P_res - P_sat) / P_res
%    - Refinement: Bisection method on ΔP → 0
%
% 2. UNDERSATURATED GOC: Bubble point equals dew point (critical point)
%    - P_bub = P_dew at GOC
%    - Detection: Saturation type remains same but ΔK approaches zero
%    - ΔK = Σ(ln K_i)² where K_i are equilibrium ratios
%    - Refinement: Secant method on ΔK → 0
%
% Algorithm Flow (from Figure 4.2):
% ---------------------------------
% 1. Start from reference depth, determine initial saturation type
% 2. Scan through depths with coarse step
% 3. At each depth:
%    a. Calculate P_res (from compositional grading)
%    b. Calculate P_bub and P_dew
%    c. Determine saturation type (Pb-dominant or Pd-dominant)
%    d. Check for GOC indicators:
%       - Saturated: Sign change in (P_res - P_sat)
%       - Undersaturated: ΔK approaching zero while P_res > P_sat
% 4. Refine GOC location:
%    - Saturated: Bisection on ΔP
%    - Undersaturated: Secant method on ΔK
%
% =========================================================================
% INPUTS:
%   comp_ref    : Reference composition [mole fractions]
%   press_ref   : Reference pressure [Pa]
%   temp        : Temperature [K] (constant for isothermal)
%   h_ref       : Reference depth [m]
%   depth_range : [h_min, h_max] - Depth range to search [m]
%   Pc, Tc, acentric, BIP, M_gmol : Fluid properties
%   vt_type     : Volume translation method (0-5)
%   vt_params   : VT parameters structure
%
% OUTPUTS:
%   GOC_depth   : Depth of GOC [m] (NaN if not found)
%   GOC_pressure: Pressure at GOC [Pa]
%   GOC_comp    : Composition at GOC [mole fractions]
%   info        : Diagnostic structure with all scan data
%
% =========================================================================

    %% Initialize
    comp_ref = comp_ref(:);
    n = length(comp_ref);
    
    % Tolerances
    tol_depth = 0.5;            % Depth tolerance [m] for convergence
    tol_DeltaK = 1e-6;          % ΔK tolerance for undersaturated GOC
    max_bisection = 30;         % Max bisection iterations
    max_secant = 30;            % Max secant iterations
    scan_step = 5;              % Coarse scan step [m]
    
    % Parse depth range
    h_min = min(depth_range);
    h_max = max(depth_range);
    
    % Determine scan direction (from reference toward search range)
    if h_ref >= h_max
        % Reference is deeper, scan upward (toward shallower)
        dh = -scan_step;
        h_start = h_ref;
        h_end = h_min;
    else
        % Reference is shallower, scan downward (toward deeper)
        dh = scan_step;
        h_start = h_ref;
        h_end = h_max;
    end
    
    depths = h_start:dh:h_end;
    n_depths = length(depths);
    
    %% Print header
    fprintf('=========================================================================\n');
    fprintf('GOC DETECTION - Nikpoor (2014) Algorithm\n');
    fprintf('=========================================================================\n');
    fprintf('Reference: h = %.1f m, P = %.2f bar, T = %.1f K\n', h_ref, press_ref/1e5, temp);
    fprintf('Scan: %.1f to %.1f m (step = %.1f m)\n', h_start, h_end, dh);
    fprintf('=========================================================================\n\n');
    
    %% Storage arrays
    info.depths = zeros(n_depths, 1);
    info.P_res = zeros(n_depths, 1);
    info.P_bub = zeros(n_depths, 1);
    info.P_dew = zeros(n_depths, 1);
    info.P_sat = zeros(n_depths, 1);
    info.Delta_P = zeros(n_depths, 1);      % (P_res - P_sat) / P_res
    info.Delta_K = zeros(n_depths, 1);      % Σ(ln K_i)²
    info.sat_type = cell(n_depths, 1);      % 'Bubble' or 'Dew'
    info.comp = zeros(n_depths, n);
    info.K_values = zeros(n_depths, n);     % K-values at each depth
    
    %% Print table header
    fprintf('%-8s %-12s %-12s %-12s %-12s %-10s %-12s %-12s\n', ...
            'Depth', 'P_res[bar]', 'P_bub[bar]', 'P_dew[bar]', 'P_sat[bar]', 'Sat Type', 'ΔP', 'ΔK');
    fprintf('%s\n', repmat('-', 1, 96));
    
    %% Scan depths
    GOC_found = false;
    GOC_type = '';              % 'saturated' or 'undersaturated'
    GOC_bracket = [];
    prev_sign = NaN;
    prev_sat_type = '';
    prev_DeltaK = NaN;
    DeltaK_decreasing_count = 0;
    
    for i = 1:n_depths
        h = depths(i);
        info.depths(i) = h;
        
        % Get composition and pressure at depth h
        try
            [comp_h, P_h, P_bub, P_dew] = main(h, h_ref, comp_ref, press_ref, temp, ...
                Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
            comp_h = comp_h(:);
        catch ME
            fprintf('Error at h = %.1f m: %s\n', h, ME.message);
            continue;
        end
        
        % Handle NaN saturation pressures
        P_bub_valid = P_bub;
        P_dew_valid = P_dew;
        if isnan(P_bub), P_bub_valid = 0; end
        if isnan(P_dew), P_dew_valid = 0; end
        
        % Determine saturation type and P_sat
        if P_bub_valid >= P_dew_valid
            sat_type = 'Bubble';
            P_sat = P_bub_valid;
        else
            sat_type = 'Dew';
            P_sat = P_dew_valid;
        end
        
        % Calculate ΔP = (P_res - P_sat) / P_res
        if P_h > 0 && P_sat > 0
            Delta_P = (P_h - P_sat) / P_h;
        else
            Delta_P = NaN;
        end
        
        % Calculate ΔK = Σ(ln K_i)² using stability analysis
        [Delta_K, K_values] = calculate_DeltaK(comp_h, P_h, temp, Pc, Tc, acentric, BIP);
        
        % Store results
        info.P_res(i) = P_h;
        info.P_bub(i) = P_bub;
        info.P_dew(i) = P_dew;
        info.P_sat(i) = P_sat;
        info.Delta_P(i) = Delta_P;
        info.Delta_K(i) = Delta_K;
        info.sat_type{i} = sat_type;
        info.comp(i, :) = comp_h';
        info.K_values(i, :) = K_values';
        
        % Print row
        fprintf('%-8.1f %-12.2f %-12.2f %-12.2f %-12.2f %-10s %-12.6f %-12.6f', ...
                h, P_h/1e5, P_bub/1e5, P_dew/1e5, P_sat/1e5, sat_type, Delta_P, Delta_K);
        
        % =====================================================================
        % GOC Detection Logic (Nikpoor 2014 Figure 4.2)
        % =====================================================================
        
        % Case 1: SATURATED GOC - Sign change in (P_res - P_sat)
        current_sign = sign(P_h - P_sat);
        
        if ~isnan(prev_sign) && current_sign ~= prev_sign && current_sign ~= 0 && prev_sign ~= 0
            % Sign change detected - SATURATED GOC
            fprintf(' <-- SATURATED GOC BRACKET');
            GOC_found = true;
            GOC_type = 'saturated';
            GOC_bracket = [depths(i-1), h];
        end
        
        % Case 2: UNDERSATURATED GOC - ΔK approaching zero while undersaturated
        % (P_res > P_sat throughout, but P_bub approaches P_dew)
        if ~GOC_found && Delta_P > 0  % Undersaturated condition
            % Check if ΔK is monotonically decreasing toward zero
            if ~isnan(prev_DeltaK)
                if Delta_K < prev_DeltaK
                    DeltaK_decreasing_count = DeltaK_decreasing_count + 1;
                else
                    DeltaK_decreasing_count = 0;
                end
                
                % Check for undersaturated GOC: ΔK approaching zero
                if Delta_K < tol_DeltaK && DeltaK_decreasing_count >= 2
                    fprintf(' <-- UNDERSATURATED GOC (ΔK → 0)');
                    GOC_found = true;
                    GOC_type = 'undersaturated';
                    GOC_bracket = [depths(i-1), h];
                end
                
                % Also check for saturation type change while undersaturated
                if ~isempty(prev_sat_type) && ~strcmp(sat_type, prev_sat_type) && Delta_P > 0.01
                    % Saturation type changed (Pb ↔ Pd dominant) while undersaturated
                    fprintf(' <-- UNDERSATURATED GOC (Type Change)');
                    GOC_found = true;
                    GOC_type = 'undersaturated';
                    GOC_bracket = [depths(i-1), h];
                end
            end
        end
        
        fprintf('\n');
        
        if GOC_found
            break;
        end
        
        % Update previous values
        prev_sign = current_sign;
        prev_sat_type = sat_type;
        prev_DeltaK = Delta_K;
    end
    
    %% Refine GOC location
    if GOC_found
        fprintf('\n=========================================================================\n');
        fprintf('GOC Type: %s\n', upper(GOC_type));
        fprintf('Bracket: [%.2f, %.2f] m\n', GOC_bracket(1), GOC_bracket(2));
        fprintf('=========================================================================\n');
        
        if strcmp(GOC_type, 'saturated')
            % Use BISECTION method for saturated GOC (ΔP → 0)
            [GOC_depth, GOC_pressure, GOC_comp] = refine_saturated_GOC(...
                GOC_bracket, h_ref, comp_ref, press_ref, temp, ...
                Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, ...
                tol_depth, max_bisection);
        else
            % Use SECANT method for undersaturated GOC (ΔK → 0)
            [GOC_depth, GOC_pressure, GOC_comp] = refine_undersaturated_GOC(...
                GOC_bracket, h_ref, comp_ref, press_ref, temp, ...
                Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, ...
                tol_depth, tol_DeltaK, max_secant);
        end
        
    else
        % No GOC found
        GOC_depth = NaN;
        GOC_pressure = NaN;
        GOC_comp = comp_ref;
        fprintf('\nNo GOC detected in range [%.1f, %.1f] m\n', h_min, h_max);
    end
    
    %% Final output
    fprintf('\n=========================================================================\n');
    fprintf('RESULT\n');
    fprintf('=========================================================================\n');
    if ~isnan(GOC_depth)
        fprintf('GOC Type:     %s\n', upper(GOC_type));
        fprintf('GOC Depth:    %.2f m\n', GOC_depth);
        fprintf('GOC Pressure: %.2f bar\n', GOC_pressure/1e5);
    else
        fprintf('No GOC found in search range.\n');
    end
    fprintf('=========================================================================\n');
    
    % Store in info
    info.GOC_depth = GOC_depth;
    info.GOC_pressure = GOC_pressure;
    info.GOC_comp = GOC_comp;
    info.GOC_type = GOC_type;
    info.GOC_bracket = GOC_bracket;
    
end


%% =========================================================================
% HELPER FUNCTION: Calculate ΔK using stability analysis
% =========================================================================
function [Delta_K, K_values] = calculate_DeltaK(comp, press, temp, Pc, Tc, acentric, BIP)
% Calculate ΔK = Σ(ln K_i)² using stability analysis to get K-values
%
% At equilibrium, K_i = y_i / x_i
% ΔK = Σ(ln K_i)² → 0 at critical point (where all K_i → 1)

    n = length(comp);
    K_values = ones(n, 1);  % Default to K=1 (critical point)
    Delta_K = 0;
    
    try
        % Use stability analysis to check if system is two-phase
        [stability, x_trial1, x_trial2, ~, ~] = stability_analysis_ssi(...
            comp, press, temp, Pc, Tc, acentric, BIP, 1e-9, 100);
        
        if stability == 2  % Unstable (two-phase)
            % Perform flash calculation to get equilibrium K-values
            % Use Rachford-Rice with trial compositions as initial guess
            try
                % Get vapor and liquid compositions from flash
                [K_flash, ~, converged] = flash_calculation_simple(...
                    comp, press, temp, Pc, Tc, acentric, BIP, x_trial1, x_trial2);
                
                if converged
                    K_values = K_flash;
                else
                    % Use Wilson K-values as fallback
                    K_values = wilsoneq(press, temp, Pc, Tc, acentric);
                end
            catch
                % Use Wilson K-values as fallback
                K_values = wilsoneq(press, temp, Pc, Tc, acentric);
            end
        else
            % Single phase - use Wilson K-values for estimation
            K_values = wilsoneq(press, temp, Pc, Tc, acentric);
        end
        
        % Calculate ΔK = Σ(ln K_i)²
        ln_K = log(K_values);
        Delta_K = sum(ln_K.^2);
        
    catch
        % If stability analysis fails, use Wilson K-values
        K_values = wilsoneq(press, temp, Pc, Tc, acentric);
        ln_K = log(K_values);
        Delta_K = sum(ln_K.^2);
    end
    
end


%% =========================================================================
% HELPER FUNCTION: Simple flash calculation for K-values
% =========================================================================
function [K, V_frac, converged] = flash_calculation_simple(z, P, T, Pc, Tc, acentric, BIP, x_init, y_init)
% Simple flash calculation to get equilibrium K-values
% Uses successive substitution

    z = z(:);
    n = length(z);
    tol = 1e-8;
    maxiter = 100;
    converged = false;
    
    % Initialize K-values from trial compositions
    x = x_init(:);
    y = y_init(:);
    
    % Ensure valid compositions
    x = max(x, 1e-15);
    y = max(y, 1e-15);
    x = x / sum(x);
    y = y / sum(y);
    
    K = y ./ x;
    K = max(K, 1e-10);
    K = min(K, 1e10);
    
    % Successive substitution
    for iter = 1:maxiter
        K_old = K;
        
        % Solve Rachford-Rice for vapor fraction
        V_frac = solve_rachford_rice(z, K);
        
        % Calculate phase compositions
        x = z ./ (1 + V_frac * (K - 1));
        y = K .* x;
        
        % Normalize
        x = x / sum(x);
        y = y / sum(y);
        
        % Calculate fugacity coefficients
        [phi_L, ~] = fugacitycoef_multicomp(x, P, T, Pc, Tc, acentric, BIP);
        [phi_V, ~] = fugacitycoef_multicomp(y, P, T, Pc, Tc, acentric, BIP);
        
        % Update K-values
        K = phi_L ./ phi_V;
        K = max(K, 1e-10);
        K = min(K, 1e10);
        
        % Check convergence
        err = max(abs(K - K_old) ./ K);
        if err < tol
            converged = true;
            break;
        end
    end
    
end


%% =========================================================================
% HELPER FUNCTION: Solve Rachford-Rice equation
% =========================================================================
function V = solve_rachford_rice(z, K)
% Solve Rachford-Rice equation for vapor fraction V
% Σ z_i(K_i - 1) / (1 + V(K_i - 1)) = 0

    z = z(:);
    K = K(:);
    
    % Bounds for V
    V_min = max(0, max((K.*z - 1)./(K - 1)));
    V_max = min(1, min((1 - z)./(1 - K)));
    
    if V_min >= V_max
        V = 0.5;
        return;
    end
    
    % Bisection
    tol = 1e-10;
    maxiter = 100;
    
    for iter = 1:maxiter
        V = (V_min + V_max) / 2;
        
        % Rachford-Rice function
        f = sum(z .* (K - 1) ./ (1 + V * (K - 1)));
        
        if abs(f) < tol
            break;
        end
        
        if f > 0
            V_min = V;
        else
            V_max = V;
        end
    end
    
end


%% =========================================================================
% HELPER FUNCTION: Refine SATURATED GOC using Bisection
% =========================================================================
function [GOC_depth, GOC_pressure, GOC_comp] = refine_saturated_GOC(...
    bracket, h_ref, comp_ref, press_ref, temp, ...
    Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, ...
    tol_depth, max_iter)
% Bisection method to find where ΔP = (P_res - P_sat) = 0

    fprintf('\nRefining SATURATED GOC using BISECTION (ΔP → 0)...\n');
    
    h_low = min(bracket);
    h_high = max(bracket);
    
    % Evaluate at boundaries
    [~, P_low, Psat_low] = eval_at_depth(h_low, h_ref, comp_ref, press_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    [~, P_high, Psat_high] = eval_at_depth(h_high, h_ref, comp_ref, press_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    
    f_low = P_low - Psat_low;
    f_high = P_high - Psat_high;
    
    fprintf('Initial: h_low=%.2f m (f=%.2e), h_high=%.2f m (f=%.2e)\n', ...
            h_low, f_low, h_high, f_high);
    
    % Bisection iterations
    for iter = 1:max_iter
        h_mid = (h_low + h_high) / 2;
        
        [comp_mid, P_mid, Psat_mid] = eval_at_depth(h_mid, h_ref, comp_ref, press_ref, temp, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
        
        f_mid = P_mid - Psat_mid;
        
        fprintf('  Iter %2d: h = %.3f m, P = %.2f bar, Psat = %.2f bar, f = %.2e\n', ...
                iter, h_mid, P_mid/1e5, Psat_mid/1e5, f_mid);
        
        % Check convergence
        if abs(h_high - h_low) < tol_depth || abs(f_mid) < 1e-3 * P_mid
            GOC_depth = h_mid;
            GOC_pressure = P_mid;
            GOC_comp = comp_mid;
            fprintf('  Converged! GOC at h = %.2f m\n', GOC_depth);
            return;
        end
        
        % Bisection update
        if f_low * f_mid < 0
            h_high = h_mid;
            f_high = f_mid;
        else
            h_low = h_mid;
            f_low = f_mid;
        end
    end
    
    % Max iterations reached
    GOC_depth = h_mid;
    GOC_pressure = P_mid;
    GOC_comp = comp_mid;
    fprintf('  Max iterations. Best estimate: h = %.2f m\n', GOC_depth);
    
end


%% =========================================================================
% HELPER FUNCTION: Refine UNDERSATURATED GOC using Secant Method
% =========================================================================
function [GOC_depth, GOC_pressure, GOC_comp] = refine_undersaturated_GOC(...
    bracket, h_ref, comp_ref, press_ref, temp, ...
    Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params, ...
    tol_depth, tol_DeltaK, max_iter)
% Secant method to find where ΔK = Σ(ln K_i)² → 0

    fprintf('\nRefining UNDERSATURATED GOC using SECANT (ΔK → 0)...\n');
    
    h_0 = bracket(1);
    h_1 = bracket(2);
    
    % Evaluate ΔK at initial points
    [comp_0, P_0, ~] = eval_at_depth(h_0, h_ref, comp_ref, press_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    [Delta_K_0, ~] = calculate_DeltaK(comp_0, P_0, temp, Pc, Tc, acentric, BIP);
    
    [comp_1, P_1, ~] = eval_at_depth(h_1, h_ref, comp_ref, press_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    [Delta_K_1, ~] = calculate_DeltaK(comp_1, P_1, temp, Pc, Tc, acentric, BIP);
    
    fprintf('Initial: h_0=%.2f m (ΔK=%.6f), h_1=%.2f m (ΔK=%.6f)\n', ...
            h_0, Delta_K_0, h_1, Delta_K_1);
    
    % Secant iterations
    for iter = 1:max_iter
        % Secant formula
        if abs(Delta_K_1 - Delta_K_0) < 1e-15
            % Denominator too small, use bisection step
            h_new = (h_0 + h_1) / 2;
        else
            h_new = h_1 - Delta_K_1 * (h_1 - h_0) / (Delta_K_1 - Delta_K_0);
        end
        
        % Ensure h_new is within reasonable bounds
        h_min = min(bracket);
        h_max = max(bracket);
        h_new = max(h_min - 10, min(h_max + 10, h_new));
        
        % Evaluate at new point
        [comp_new, P_new, ~] = eval_at_depth(h_new, h_ref, comp_ref, press_ref, temp, ...
            Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
        [Delta_K_new, ~] = calculate_DeltaK(comp_new, P_new, temp, Pc, Tc, acentric, BIP);
        
        fprintf('  Iter %2d: h = %.3f m, P = %.2f bar, ΔK = %.6f\n', ...
                iter, h_new, P_new/1e5, Delta_K_new);
        
        % Check convergence
        if abs(h_new - h_1) < tol_depth || Delta_K_new < tol_DeltaK
            GOC_depth = h_new;
            GOC_pressure = P_new;
            GOC_comp = comp_new;
            fprintf('  Converged! GOC at h = %.2f m\n', GOC_depth);
            return;
        end
        
        % Update for next iteration
        h_0 = h_1;
        Delta_K_0 = Delta_K_1;
        h_1 = h_new;
        Delta_K_1 = Delta_K_new;
    end
    
    % Max iterations reached
    GOC_depth = h_new;
    GOC_pressure = P_new;
    GOC_comp = comp_new;
    fprintf('  Max iterations. Best estimate: h = %.2f m\n', GOC_depth);
    
end


%% =========================================================================
% HELPER FUNCTION: Evaluate composition, pressure, and Psat at depth
% =========================================================================
function [comp_h, P_h, Psat] = eval_at_depth(h, h_ref, comp_ref, press_ref, temp, ...
    Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params)
% Evaluate composition, pressure, and saturation pressure at depth h

    [comp_h, P_h, P_bub, P_dew] = main(h, h_ref, comp_ref, press_ref, temp, ...
        Pc, Tc, acentric, BIP, M_gmol, vt_type, vt_params);
    comp_h = comp_h(:);
    
    % Handle NaN
    if isnan(P_bub), P_bub = 0; end
    if isnan(P_dew), P_dew = 0; end
    
    % Saturation pressure
    Psat = max(P_bub, P_dew);
    
end