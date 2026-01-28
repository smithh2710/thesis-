function [comp_h, press_h, temp_h] = main_hasse(h_target, h_ref, comp_ref, P_ref, T_ref, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref)
% Step-wise integration from h_ref to h_target

    step_size = 1;  % 1 meter steps
    
    if h_target > h_ref
        depths = h_ref:step_size:h_target;
    else
        depths = h_ref:-step_size:h_target;
    end
    
    % Initialize
    comp_current = comp_ref;
    P_current = P_ref;
    T_current = T_ref;
    
    for i = 2:length(depths)
        h_prev = depths(i-1);
        h_curr = depths(i);
        
        % Single step calculation
        [comp_current, P_current, T_current] = main_hasse_single_step(...
            h_curr, h_prev, comp_current, P_current, T_current, dTdh, ...
            Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);
    end
    
    comp_h = comp_current;
    press_h = P_current;
    temp_h = T_current;
    
    fprintf('Step-wise result at %.0f m: P = %.2f bar\n', h_target, press_h/1e5);

end


function [comp_h, press_h, temp_h] = main_hasse_single_step(h, h_ref, comp_ref, P, T, dTdh, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref)
% Single step Haase calculation (small delta_h)

    R = 8.3144598;
    g = 9.80665;
    n = length(comp_ref);

    comp_ref = comp_ref(:);
    comp_ref = comp_ref / sum(comp_ref);
    M_gmol = M_gmol(:);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);

    delta_h = h - h_ref;
    delta_T = dTdh * delta_h;
    temp_h = T + delta_T;

    % Reference fugacities
    [fugcoef_ref, ~] = fugacitycoef_multicomp(comp_ref, P, T, Pc, Tc, acentric, BIP);
    f_ref = fugcoef_ref .* comp_ref * P;

    % Enthalpy at current step
    [H_mix_specific, H_partial_specific, ~, ~, ~] = calculate_absolute_enthalpy(...
        temp_h, P, comp_ref, Pc, Tc, acentric, BIP, M_gmol, Cp_coeffs, H_ig_ref);

    % Compute terms
    term2 = zeros(n, 1);
    term3 = zeros(n, 1);

    for i = 1:n
        M_i_kg = M_gmol(i) / 1000;
        term2(i) = M_i_kg * g * delta_h / (R * T);
        
        enthalpy_diff = H_mix_specific - H_partial_specific(i);
        term3(i) = M_gmol(i) * enthalpy_diff * delta_T / (R * T^2);
    end

    % Setup parameters
    params = struct();
    params.f_ref = f_ref;
    params.term2 = term2;
    params.term3 = term3;
    params.T_ref = T;
    params.Pc = Pc;
    params.Tc = Tc;
    params.acentric = acentric;
    params.BIP = BIP;
    params.n = n;

    initial_guess = [comp_ref; P];

    options = optimoptions('fsolve', ...
        'Display', 'none', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'MaxIterations', 200);

    [solution, ~, ~] = fsolve(@(x) residual_haase(x, params), initial_guess, options);

    comp_h = solution(1:n);
    press_h = solution(end);

    comp_h = max(comp_h, 1e-15);
    comp_h = comp_h / sum(comp_h);

end


function F = residual_haase(x, params)

    f_ref = params.f_ref;
    term2 = params.term2;
    term3 = params.term3;
    T_ref = params.T_ref;
    Pc = params.Pc;
    Tc = params.Tc;
    acentric = params.acentric;
    BIP = params.BIP;
    n = params.n;

    z_h = x(1:n);
    P_h = x(end);

    z_h = max(z_h, 1e-15);
    z_h_norm = z_h / sum(z_h);
    P_h = max(P_h, 1e2);

    [fugcoef_h, ~] = fugacitycoef_multicomp(z_h_norm, P_h, T_ref, Pc, Tc, acentric, BIP);

    F = zeros(n+1, 1);

    for i = 1:n
        f_i_h = fugcoef_h(i) * z_h_norm(i) * P_h;
        F(i) = log(f_i_h / f_ref(i)) - term2(i) + term3(i);
    end

    F(n+1) = sum(z_h) - 1;

end
