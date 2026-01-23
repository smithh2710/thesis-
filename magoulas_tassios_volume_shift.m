function [s, c, mix] = magoulas_tassios_volume_shift(T, Pc, Tc, acentric, comp, Vc)
% MAGOULAS_TASSIOS_VOLUME_SHIFT - Volume translation for PR-EOS using Magoulas & Tassios (1990)

% Reference:
%   Magoulas, K. and Tassios, D. (1990). "Thermophysical properties of 
%   n-Alkanes from C1 to C20 and their prediction"
%   Fluid Phase Equilibria, Vol. 56, pp. 119-140.

% Equations:
%   c = -[c_m0 + (δc - c_m0) × exp(β × |1 - T/Tc|)]           
%   c_m0 = (R×Tc/Pc) × (k_m0 + k_m1×ω + k_m2×ω² + k_m3×ω³ + k_m4×ω⁴)  
%   β = l_m0 + l_m1×ω²                                         
%   δc = (R×Tc/Pc) × (Zc_PR - Zc_exp)                        
%   Zc_exp = 0.289 - 0.0701×ω - 0.0207×ω²                    
%
% Model Parameters:
%   k_m0 = -0.014471,  k_m1 = 0.067498,  k_m2 = -0.084852
%   k_m3 = 0.067298,   k_m4 = -0.017366
%   l_m0 = -10.2447,   l_m1 = -28.6312
%   Zc_PR = 0.3074
%
% -------------------------------------------------------------------------
% Inputs:
%   T        : Temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   comp     : (Optional) Mole fractions [-] for mixture calculation
%   Vc       : (Optional) Critical volumes [cm³/mol] for seta mixing rule
%              If not provided, estimated from Zc correlation
%
% Outputs:
%   s   : Dimensionless volume shift parameters [-] (s = c/b)
%   c   : Volume shift parameters [m³/mol]
%   mix : Structure with mixture and component properties:
%         .c_mix      - Mixture volume translation [m³/mol]
%         .c_mix_cm3  - Mixture volume translation [cm³/mol]
%         .s_mix      - Mixture dimensionless shift [-]
%         .seta       - Surface fractions [-]
%         .c_cm3      - Individual c values [cm³/mol]
%         .c_m0       - Individual c_m0 values [m³/mol]
%         .delta_c    - Individual δc values [m³/mol]
%         .beta       - Individual β values [-]
%         .Zc_exp     - Estimated experimental Zc [-]


    % Constants
    R = 8.3144598;           % J/(mol·K)
    Omega_b = 0.07780;       % PR-EOS b parameter constant
    Zc_PR = 0.3074;          % PR-EOS critical compressibility
    
    % Model parameters (from Shi 2017 thesis, Table 2.1)
    k_m0 = -0.014471;
    k_m1 =  0.067498;
    k_m2 = -0.084852;
    k_m3 =  0.067298;
    k_m4 = -0.017366;
    
    l_m0 = -10.2447;
    l_m1 = -28.6312;
    
    % Ensure column vectors
    nc = length(Pc);
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    
    % Handle optional composition input
    if nargin >= 5 && ~isempty(comp)
        comp = comp(:);
    else
        comp = ones(nc, 1) / nc;  % Default equal molar
    end
    
    % Handle optional Vc input
    if nargin >= 6 && ~isempty(Vc)
        Vc = Vc(:);
    else
        % Estimate Vc from Zc_exp correlation (Eq. 21)
        Zc_exp_est = 0.289 - 0.0701 * acentric - 0.0207 * acentric.^2;
        Vc = Zc_exp_est .* R .* Tc ./ Pc * 1e6;  % [cm³/mol]
    end
    
    % Preallocate
    c = zeros(nc, 1);
    c_cm3 = zeros(nc, 1);
    s = zeros(nc, 1);
    b = zeros(nc, 1);
    c_m0 = zeros(nc, 1);
    delta_c = zeros(nc, 1);
    beta = zeros(nc, 1);
    Zc_exp = zeros(nc, 1);
    
    % Calculate volume translation for each component
    for i = 1:nc
        w = acentric(i);
        Tc_i = Tc(i);
        Pc_i = Pc(i);
        
        % Reduced temperature term
        Tr = T / Tc_i;
        tau = abs(1 - Tr);  % |1 - T/Tc|
        
        % Eq. 21: Experimental critical compressibility
        Zc_exp(i) = 0.289 - 0.0701*w - 0.0207*w^2;
        
        % Eq. 20: Critical point correction δc
        delta_c(i) = (R * Tc_i / Pc_i) * (Zc_PR - Zc_exp(i));  % [m³/mol]
        
        % Eq. 18: c_m0 (volume shift contribution at critical)
        k_poly = k_m0 + k_m1*w + k_m2*w^2 + k_m3*w^3 + k_m4*w^4;
        c_m0(i) = (R * Tc_i / Pc_i) * k_poly;  % [m³/mol]
        
        % Eq. 19: β (exponential decay parameter)
        beta(i) = l_m0 + l_m1*w^2;
        
        % Eq. 17: Volume translation c
        % c = -[c_m0 + (δc - c_m0) × exp(β × |1 - T/Tc|)]
        c(i) = (c_m0(i) + (delta_c(i) - c_m0(i)) * exp(beta(i) * tau));  % [m³/mol]
        
        % Convert to cm³/mol
        c_cm3(i) = c(i) * 1e6;
        
        % Calculate b parameter for dimensionless shift
        b(i) = Omega_b * R * Tc_i / Pc_i;  % [m³/mol]
        
        % Dimensionless shift
        s(i) = c(i) / b(i);
    end
    
    % Calculate mixture properties
    mix = struct();
    mix.c_cm3 = c_cm3;
    mix.c_m0 = c_m0;
    mix.delta_c = delta_c;
    mix.beta = beta;
    mix.Zc_exp = Zc_exp;
    
    % Surface fraction (seta) mixing rule - same as Peneloux/Baled
    % seta_i = (z_i × Vc_i^(2/3)) / sum(z_j × Vc_j^(2/3))
    Vc_23 = Vc.^(2/3);
    sum_zVc23 = sum(comp .* Vc_23);
    
    if sum_zVc23 > 0
        seta = (comp .* Vc_23) / sum_zVc23;
    else
        seta = comp;  % Fallback to molar average
    end
    
    % Mixture volume translation (surface fraction weighted)
    c_mix = sum(seta .* c);              % [m³/mol]
    c_mix_cm3 = c_mix * 1e6;             % [cm³/mol]
    
    % Mixture b parameter (molar average)
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
    mix.Zc_PR = Zc_PR;                   % [-]
    
end