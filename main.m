function [comp_h, press_h, pressbub_h, pressdew_h] = main(h, h_ref, comp_ref, press_ref, temp, pressc, tempc, acentric, BIP, M_gmol, vt_method, vt_params)
% MAIN - Isothermal Compositional Grading with Volume Translation
% =========================================================================
% Solves the isothermal compositional grading equations (Schulte 1980):
%   ln(φ_i^h z_i^h P^h) - ln(φ_i^h° z_i^h° P^h°) = M_i g (h - h°) / (RT)
%   Σ z_i = 1
%
% INPUTS:
%   h         : Target depth [m]
%   h_ref     : Reference depth [m]
%   comp_ref  : Reference composition [mole fractions]
%   press_ref : Reference pressure [Pa]
%   temp      : Temperature [K] (constant for isothermal)
%   pressc    : Critical pressures [Pa]
%   tempc     : Critical temperatures [K]
%   acentric  : Acentric factors [-]
%   BIP       : Binary interaction parameter matrix
%   M_gmol    : Molecular weights [g/mol]
%   vt_method : Volume translation method:
%               0  - None
%               1  - Peneloux (1982)
%               2  - Magoulas-Tassios (1990)
%               3  - Ungerer-Batut (1997)
%               4  - Baled (2012)
%               5  - Abudour (2012)
%               6  - CUSTOM: User-provided c values in vt_params.c_custom
%               OR: Direct vector of c values [m³/mol] (auto-detected)
%   vt_params : Structure with additional parameters:
%               .Vc         - Critical volumes [cm³/mol] (for some methods)
%               .components - Cell array of component names (for Baled/Abudour)
%               .c_custom   - Custom volume translation [cm³/mol] (for method 6)
%
% OUTPUTS:
%   comp_h    : Composition at depth h [mole fractions]
%   press_h   : Pressure at depth h [Pa]
%   pressbub_h: Bubble point pressure at depth h [Pa]
%   pressdew_h: Dew point pressure at depth h [Pa]
%
% =========================================================================
% VOLUME TRANSLATION METHOD 6 - CUSTOM INPUT:
%   Set vt_method = 6 and provide vt_params.c_custom as a vector of 
%   volume corrections in [cm³/mol] for each component.
%
%   Example (Reservoir 1 from Pedersen et al. 2015, Table 7):
%     vt_params.c_custom = [-4.23; -1.91; -5.20; -5.79; -6.35; ...];  % cm³/mol
%     [comp_h, P_h, Pb, Pd] = main(h, h_ref, comp, P_ref, T, Pc, Tc, w, BIP, M, 6, vt_params);
%
% =========================================================================

    % Constants
    R = 8.3144598;        % J/(mol·K)
    g = 9.80665;          % m/s²
    M = M_gmol / 1000;    % Convert to kg/mol

    tol = 1e-10;
    maxiter = 1500;
    n = length(comp_ref);
    
    % Ensure column vectors
    comp_ref = comp_ref(:);
    M = M(:);
    M_gmol = M_gmol(:);
    pressc = pressc(:);
    tempc = tempc(:);
    acentric = acentric(:);

    % Handle optional inputs
    if nargin < 12
        vt_params = struct();
    end
    if nargin < 11
        vt_method = 0;
    end

    % Parse volume translation method
    [vt_type, vt_opts] = parse_vt_input(vt_method, vt_params, n, temp, pressc, tempc, acentric, M_gmol);
    
    % Calculate volume translation at reference conditions
    c_ref = calc_volume_translation(comp_ref, press_ref, temp, pressc, tempc, acentric, M_gmol, vt_type, vt_opts);

    % Calculate fugacity coefficients at reference
    [fugcoef_ref_eos, ~] = fugacitycoef_multicomp(comp_ref, press_ref, temp, pressc, tempc, acentric, BIP);
    
    % Apply volume translation correction at reference
    if vt_type > 0 || vt_type == -1  % Include custom VT (type -1)
        fugcoef_ref = fugcoef_ref_eos .* exp(-c_ref * press_ref / (R * temp));
    else
        fugcoef_ref = fugcoef_ref_eos;
    end
    
    % Reference fugacity
    f_ref = fugcoef_ref .* comp_ref * press_ref;

    % Fugacity at depth h (Eq. 3.27)
    f_h = f_ref .* exp(( M * g * (h - h_ref)) / (R * temp));

    % Initial guess
    initial_guess = [comp_ref; press_ref];

    % Solve for composition and pressure at depth h
    fun = @(x) residual_fugacity(x(1:n), x(end), f_h, temp, pressc, tempc, acentric, BIP, M_gmol, vt_type, vt_opts, R);
    
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-15, 'StepTolerance', 1e-15);
    solution = fsolve(fun, initial_guess, options);

    comp_h = solution(1:n);
    press_h = solution(end);

    % Bubble pressure at h
    try 
        pressbub_ini_h = 270e5;
        [pressbub_h, ~] = pressbub_multicomp_newton(comp_h, pressbub_ini_h, temp, pressc, tempc, acentric, BIP, tol, maxiter);
        if ~isreal(pressbub_h) || pressbub_h <= 0 || ~isfinite(pressbub_h)
           pressbub_h = NaN;
        end
    catch 
       pressbub_h = NaN; 
    end   

    % Dew pressure at h
    try 
        pressdew_ini_h = 270e5;
        [pressdew_h, ~] = pressdew_multicomp_newton(comp_h, pressdew_ini_h, temp, pressc, tempc, acentric, BIP, tol, maxiter);
        if ~isreal(pressdew_h) || pressdew_h <= 0 || ~isfinite(pressdew_h)
           pressdew_h = NaN;
        end
    catch 
       pressdew_h = NaN; 
    end 

end


function F = residual_fugacity(comp_h, press_h, f_h, temp, pressc, tempc, acentric, BIP, M_gmol, vt_type, vt_opts, R)
% RESIDUAL_FUGACITY - Residual function for compositional grading
%
% Volume translation is calculated at each iteration to ensure consistency,
% especially for pressure/density-dependent methods like Abudour.

    n = length(f_h);
    
    % Get EOS fugacity coefficients
    [fugcoef_h_eos, ~] = fugacitycoef_multicomp(comp_h, press_h, temp, pressc, tempc, acentric, BIP);
    
    % Calculate volume translation at current conditions
    if vt_type > 0 || vt_type == -1  % Include custom VT (type -1)
        c_shift = calc_volume_translation(comp_h, press_h, temp, pressc, tempc, acentric, M_gmol, vt_type, vt_opts);
        fugcoef_h = fugcoef_h_eos .* exp(-c_shift * press_h / (R * temp));
    else
        fugcoef_h = fugcoef_h_eos;
    end
    
    % Residual equations
    F = zeros(n+1, 1);
    for i = 1:n
        F(i) = comp_h(i) * press_h * fugcoef_h(i) - f_h(i);
    end
    F(n+1) = sum(comp_h) - 1;
end


function c = calc_volume_translation(comp, press, temp, Pc, Tc, acentric, M_gmol, vt_type, vt_opts)
% CALC_VOLUME_TRANSLATION - Calculate volume translation for any method
%
% Inputs:
%   comp     : Composition [mole fractions]
%   press    : Pressure [Pa]
%   temp     : Temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   M_gmol   : Molecular weights [g/mol]
%   vt_type  : VT method number (0-6 or -1 for direct input)
%   vt_opts  : Structure with additional parameters
%
% Output:
%   c        : Volume translation [m³/mol] for each component

    n = length(comp);
    
    switch vt_type
        case 0  % No volume translation
            c = zeros(n, 1);
            
        case -1  % Direct c input (auto-detected vector or c_custom in cm³/mol)
            c = vt_opts.c_direct;  % Already in m³/mol
            
        case 1  % Peneloux (T-independent)
            [~, c] = peneloux_volume_shift(Pc, Tc, acentric);
            
        case 2  % Magoulas-Tassios (T-dependent)
            [~, c] = magoulas_tassios_volume_shift(temp, Pc, Tc, acentric);
            
        case 3  % Ungerer-Batut (T,MW-dependent)
            [~, c] = ungerer_batut_volume_shift(temp, Pc, Tc, acentric, M_gmol);
            
        case 4  % Baled (T-dependent)
            components = vt_opts.components;
            [~, c] = baled_volume_shift(temp, Pc, Tc, acentric, M_gmol, components);
            
        case 5  % Abudour (T,P,ρ-dependent)
            comp = comp(:);
            Vc = vt_opts.Vc;
            components = vt_opts.components;
            BIP_local = zeros(n);  % Use zero BIP for internal VT calculation
            [~, c] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, Vc, components, BIP_local, false);
            
        case 6  % CUSTOM: User-provided c values
            c = vt_opts.c_direct;  % Already converted to m³/mol in parse_vt_input
            
        otherwise
            c = zeros(n, 1);
    end
    
    c = c(:);  % Ensure column vector
end


function [vt_type, vt_opts] = parse_vt_input(vt_method, vt_params, n, temp, Pc, Tc, acentric, M_gmol)
% PARSE_VT_INPUT - Parse volume translation method and options
%
% VT Methods for PR-EOS:
%   0  - None
%   1  - Peneloux (1982) - adapted for PR
%   2  - Magoulas-Tassios (1990)
%   3  - Ungerer-Batut (1997)
%   4  - Baled (2012)
%   5  - Abudour (2012)
%   6  - CUSTOM: User-provided c values in vt_params.c_custom [cm³/mol]
%   -1 - Direct c input (auto-detected when vt_method is a vector)

    vt_opts = struct();
    R = 8.3144598;
    
    % Extract optional parameters with defaults
    if isfield(vt_params, 'Vc')
        vt_opts.Vc = vt_params.Vc(:);
    else
        % Estimate Vc from Zc correlation
        Zc_est = 0.2905 - 0.085 * acentric;
        vt_opts.Vc = Zc_est .* R .* Tc ./ Pc * 1e6;  % cm³/mol
    end
    
    if isfield(vt_params, 'components')
        vt_opts.components = vt_params.components;
    else
        vt_opts.components = cell(n, 1);
    end
    
    % Store parameters needed for VT calculations
    vt_opts.Pc = Pc;
    vt_opts.Tc = Tc;
    vt_opts.acentric = acentric;
    vt_opts.M_gmol = M_gmol;
    vt_opts.temp = temp;
    vt_opts.n = n;
    
    % Handle empty input
    if isempty(vt_method)
        vt_type = 0;
        return;
    end
    
    % =====================================================================
    % CHECK FOR DIRECT VECTOR INPUT (auto-detect)
    % If vt_method is a numeric vector of length n, treat as direct c input
    % =====================================================================
    if isnumeric(vt_method) && length(vt_method) == n
        vt_type = -1;
        % Assume input is in cm³/mol, convert to m³/mol
        vt_opts.c_direct = vt_method(:) * 1e-6;
        return;
    end
    
    % Convert string to number if needed
    if ischar(vt_method) || isstring(vt_method)
        vt_type = get_vt_method_number(vt_method);
    else
        vt_type = vt_method;
    end
    
    % =====================================================================
    % HANDLE METHOD 6: CUSTOM VT FROM vt_params.c_custom
    % =====================================================================
    if vt_type == 6
        if isfield(vt_params, 'c_custom') && ~isempty(vt_params.c_custom)
            c_custom = vt_params.c_custom(:);
            if length(c_custom) ~= n
                error('c_custom must have %d elements (one per component)', n);
            end
            % Input is in cm³/mol, convert to m³/mol
            vt_opts.c_direct = c_custom * 1e-6;
        else
            error('For vt_method = 6, you must provide vt_params.c_custom [cm³/mol]');
        end
        return;
    end
    
    % Validate
    if vt_type < 0 || vt_type > 6
        warning('Unknown VT method %d. No volume translation applied.', vt_type);
        vt_type = 0;
    end
end


function num = get_vt_method_number(name)
% GET_VT_METHOD_NUMBER - Convert string method name to number
%
% PR-EOS Volume Translation Methods:
%   0 - None
%   1 - Peneloux (adapted for PR)
%   2 - Magoulas-Tassios
%   3 - Ungerer-Batut
%   4 - Baled
%   5 - Abudour
%   6 - Custom

    name = lower(char(name));
    
    switch name
        case {'none', 'no', '0', ''}
            num = 0;
        case {'peneloux', 'pen', '1'}
            num = 1;
        case {'magoulas_tassios', 'magoulas', 'tassios', 'mt', '2'}
            num = 2;
        case {'ungerer_batut', 'ungerer', 'batut', 'ub', '3'}
            num = 3;
        case {'baled', 'bal', '4'}
            num = 4;
        case {'abudour', 'abu', '5'}
            num = 5;
        case {'custom', 'user', '6'}
            num = 6;
        otherwise
            warning('Unknown VT method: %s. No volume translation applied.', name);
            num = 0;
    end
end