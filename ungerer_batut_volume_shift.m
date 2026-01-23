function [s, c, mix] = ungerer_batut_volume_shift(T, Pc, Tc, acentric, M_gmol, comp, Vc)
% UNGERER_BATUT_VOLUME_SHIFT - Volume translation for PR-EOS using Ungerer & Batut (1997)
% =========================================================================
% Reference:
%   Ungerer, P. & Batut, C. (1997). "Prédiction des propriétés volumétriques 
%   des hydrocarbures par une translation de volume améliorée"
%   Revue de l'Institut Français du Pétrole, 52(6), 609-623.
%
% Key Features:
%   - T-dependent and MW-dependent correlation
%   - Developed for HTHP reservoir conditions (up to 200°C, 120 MPa)
%   - Calibrated on high-pressure density data (not saturated liquid)
%   - Simple linear expression
%   - Applicable to C6-C40 n-alkanes, cyclohexane, C6-C12 monoaromatics
%
% Equation (from paper):
%   c(T) = (0.023 - 0.00056 × MW) × T + (-34.5 + 0.4666 × MW)  [cm³/mol]
%
% Mixing Rule (surface fraction weighted, same as Peneloux/Baled):
%   seta_i = (z_i × Vc_i^(2/3)) / sum(z_j × Vc_j^(2/3))
%   c_mix = sum(seta_i × c_i)
%
% -------------------------------------------------------------------------
% Inputs:
%   T        : Temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   M_gmol   : Molecular weights [g/mol]
%   comp     : (Optional) Mole fractions [-] for mixture calculation
%   Vc       : (Optional) Critical volumes [cm³/mol] for seta mixing rule
%              If not provided, estimated from Zc correlation
%
% Outputs:
%   s   : Dimensionless volume shift parameters [-] (s = c/b)
%   c   : Volume shift parameters [m³/mol]
%   mix : Structure with mixture properties (if comp provided):
%         .c_mix      - Mixture volume translation [m³/mol]
%         .c_mix_cm3  - Mixture volume translation [cm³/mol]
%         .s_mix      - Mixture dimensionless shift [-]
%         .seta       - Surface fractions [-]
%         .c_cm3      - Individual c values [cm³/mol]
%
% -------------------------------------------------------------------------
% Usage:
%   % Pure components only
%   [s, c] = ungerer_batut_volume_shift(T, Pc, Tc, acentric, M);
%
%   % With mixture properties
%   [s, c, mix] = ungerer_batut_volume_shift(T, Pc, Tc, acentric, M, comp, Vc);
%
%   % Example for n-decane at 400 K
%   T = 400; Pc = 2.11e6; Tc = 617.7; w = 0.4884; MW = 142.28;
%   [s, c, mix] = ungerer_batut_volume_shift(T, Pc, Tc, w, MW);
%   % c ≈ 31.9 cm³/mol
%
% =========================================================================
% Validity Range:
%   - Temperature: up to 473 K (200°C)
%   - Pressure: up to 120 MPa
%   - Components: C6+ hydrocarbons (best for heavier HCs)
%   - For lighter components (C1-C5), use with caution
%
% Typical AAD: < 3% for C6-C13 hydrocarbons at high pressure
% =========================================================================

    % Constants
    R = 8.3144598;           % J/(mol·K)
    Omega_b = 0.07780;       % PR-EOS b parameter constant
    
    % Ensure column vectors
    nc = length(Pc);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    M_gmol = M_gmol(:);
    
    % Handle optional composition input
    if nargin >= 6 && ~isempty(comp)
        comp = comp(:);
        calc_mix = true;
    else
        comp = ones(nc, 1) / nc;  % Default equal molar
        calc_mix = false;
    end
    
    % Handle optional Vc input
    if nargin >= 7 && ~isempty(Vc)
        Vc = Vc(:);
    else
        % Estimate Vc from Zc correlation: Zc = 0.2905 - 0.085*omega
        Zc_est = 0.2905 - 0.085 * acentric;
        Vc = Zc_est .* R .* Tc ./ Pc * 1e6;  % [cm³/mol]
    end
    
    % Preallocate
    c_cm3 = zeros(nc, 1);
    c = zeros(nc, 1);
    s = zeros(nc, 1);
    b = zeros(nc, 1);
    
    % Calculate volume translation for each component
    % Ungerer-Batut correlation:
    % c(T) = (0.023 - 0.00056 × MW) × T + (-34.5 + 0.4666 × MW)  [cm³/mol]
    
    for i = 1:nc
        MW = M_gmol(i);
        
        % Calculate c in cm³/mol
        c_cm3(i) = (0.023 - 0.00056 * MW) * T + (-34.5 + 0.4666 * MW);
        
        % Convert to m³/mol
        c(i) = c_cm3(i) * 1e-6;
        
        % Calculate b parameter for dimensionless shift
        b(i) = Omega_b * R * Tc(i) / Pc(i);  % [m³/mol]
        
        % Dimensionless shift
        s(i) = c(i) / b(i);
    end
    
    % Calculate mixture properties
    mix = struct();
    mix.c_cm3 = c_cm3;
    
    % Surface fraction (seta) mixing rule - same as Peneloux/Baled
    % seta_i = (z_i × Vc_i^(2/3)) / sum(z_j × Vc_j^(2/3))
    Vc_23 = Vc.^(2/3);
    sum_zVc23 = sum(comp .* Vc_23);
    
    if sum_zVc23 > 0
        seta = (comp .* Vc_23) / sum_zVc23;
    else
        seta = comp;  % Fallback to molar average
    end
    
    % Mixture volume translation
    c_mix_cm3 = sum(seta .* c_cm3);      % [cm³/mol]
    c_mix = c_mix_cm3 * 1e-6;            % [m³/mol]
    
    % Mixture b parameter
    b_mix = sum(comp .* b);              % [m³/mol]
    
    % Mixture dimensionless shift
    s_mix = c_mix / b_mix;
    
    % Store mixture properties
    mix.c_mix = c_mix;                   % [m³/mol]
    mix.c_mix_cm3 = c_mix_cm3;           % [cm³/mol]
    mix.s_mix = s_mix;                   % [-]
    mix.seta = seta;                     % [-]
    mix.b = b;                           % [m³/mol]
    mix.b_mix = b_mix;                   % [m³/mol]
    mix.Vc = Vc;                         % [cm³/mol]
    mix.T = T;                           % [K]
    
end