function [comp_h, press_h, temp_h, pressbub_h, pressdew_h] = main_nonisothermal(h, h_ref, comp_ref, press_ref, temp_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref)
% ... [keep existing header comments]

    R = 8.3144598;
    g = 9.80665;

    n = length(comp_ref);

    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);
    M_gmol = M_gmol(:);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    H_ig_ref = H_ig_ref(:);

    M_kgmol = M_gmol / 1000;

    delta_h = h - h_ref;
    temp_h = temp_ref + dTdh * delta_h;

    tol = 1e-10;
    maxiter = 1500;

    % Calculate fugacity at reference conditions
    [fugcoef_ref, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp_ref, Pc, Tc, acentric, BIP);
    f_ref = fugcoef_ref .* comp_ref * press_ref;

    % Initial guess
    initial_guess = [comp_ref; press_ref];

    % Create parameter structure
    params.f_ref = f_ref;
    params.temp_ref = temp_ref;
    params.temp_h = temp_h;
    params.delta_h = delta_h;
    params.Pc = Pc;
    params.Tc = Tc;
    params.acentric = acentric;
    params.BIP = BIP;
    params.M_gmol = M_gmol;
    params.M_kgmol = M_kgmol;
    params.Cp_coeffs = Cp_coeffs;
    params.H_ig_ref = H_ig_ref;
    params.R = R;
    params.g = g;
    params.n = n;
    params.comp_ref = comp_ref;       % ADD THIS
    params.press_ref = press_ref;     % ADD THIS

    % Solve
    residual_fun = @(x) residual_haase_corrected(x, params);

    options = optimoptions('fsolve', ...
        'Display', 'none', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12, ...
        'MaxIterations', maxiter, ...
        'MaxFunctionEvaluations', maxiter * (n+1));

    [solution, ~, exitflag] = fsolve(residual_fun, initial_guess, options);

    if exitflag <= 0
        warning('main_nonisothermal: fsolve did not converge (exitflag = %d)', exitflag);
    end

    comp_h = solution(1:n);
    press_h = solution(end);

    comp_h = max(comp_h, 0);
    comp_h = comp_h / sum(comp_h);

    % Bubble point
    try
        pressbub_ini = 260e5;
        [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, pressbub_ini, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        if ~isreal(pressbub_h) || pressbub_h <= 0 || ~isfinite(pressbub_h)
            pressbub_h = NaN;
        end
    catch
        pressbub_h = NaN;
    end

    % Dew point
    try
        pressdew_ini = 260e5;
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, pressdew_ini, temp_h, Pc, Tc, acentric, BIP, tol, maxiter);
        if ~isreal(pressdew_h) || pressdew_h <= 0 || ~isfinite(pressdew_h)
            pressdew_h = NaN;
        end
    catch
        pressdew_h = NaN;
    end

end


function F = residual_haase_corrected(x, params)
% CORRECTED Haase residual function
%
% Key fix: Enthalpy terms calculated at REFERENCE conditions, not target conditions
% This matches the Pedersen et al. (2015) formulation more closely

    % Extract parameters
    f_ref = params.f_ref;
    temp_ref = params.temp_ref;
    temp_h = params.temp_h;
    delta_h = params.delta_h;
    Pc = params.Pc;
    Tc = params.Tc;
    acentric = params.acentric;
    BIP = params.BIP;
    M_gmol = params.M_gmol;
    M_kgmol = params.M_kgmol;
    Cp_coeffs = params.Cp_coeffs;
    H_ig_ref = params.H_ig_ref;
    R = params.R;
    g = params.g;
    n = params.n;
    comp_ref = params.comp_ref;
    press_ref = params.press_ref;

    % Extract current trial values
    z_h = x(1:n);
    P_h = x(end);

    z_h = max(z_h, 1e-15);
    z_h_norm = z_h / sum(z_h);
    P_h = max(P_h, 1e5);

    %% Calculate enthalpies at REFERENCE conditions (not target!)
    % This is a key point - Pedersen uses reference state enthalpies
    [H_abs_ref, H_mix_ref, ~, ~, ~] = calculate_absolute_enthalpy(...
        temp_ref, press_ref, comp_ref, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);

    % Mixture molecular weight at reference [g/mol]
    M_mix_ref = sum(comp_ref .* M_gmol);
    
    % Specific enthalpies at reference [J/g]
    H_specific_mix_ref = H_mix_ref / M_mix_ref;
    H_specific_i_ref = H_abs_ref ./ M_gmol;

    %% Calculate fugacity coefficients at TARGET conditions
    [fugcoef_h, ~] = fugacitycoef_multicomp(z_h_norm, P_h, temp_h, Pc, Tc, acentric, BIP);
    f_h_current = fugcoef_h .* z_h_norm * P_h;

    %% Calculate target fugacity using Haase equation
    delta_T = temp_h - temp_ref;
    f_h_target = zeros(n, 1);

    for i = 1:n
        % Gravitational term (POSITIVE for deeper locations)
        grav_term = (M_kgmol(i) * g * delta_h) / (R * temp_h);

        % Thermal term using REFERENCE enthalpies
        % (H/M - H_i/M_i) at reference conditions
        enthalpy_diff = H_specific_mix_ref - H_specific_i_ref(i);  % [J/g]
        
        % Haase thermal term: M_i * (H/M - H_i/M_i) * dT / (R * T * T_ref)
        % Units: [g/mol] * [J/g] * [K] / ([J/(mol·K)] * [K] * [K]) = dimensionless
        thermal_term = (M_gmol(i) * enthalpy_diff * delta_T) / (R * temp_h * temp_ref);

        % Target fugacity: add gravity, subtract thermal
        ln_f_target = log(f_ref(i)) + grav_term - thermal_term;
        f_h_target(i) = exp(ln_f_target);
    end

    %% Residual equations
    F = zeros(n+1, 1);
    F(1:n) = f_h_current - f_h_target;
    F(n+1) = sum(z_h) - 1;

end