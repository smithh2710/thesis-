function [H_abs, H_mix, H_ig_specific, H_res_specific, H_abs_specific] = calculate_absolute_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, Mw, Cp_coeffs, H_ig_ref)

%   H_i^abs = H_i^ig + H_i^res                              (Eq. 7)
%   H^abs = sum(z_i * H_i^abs)                              (Eq. 8)
%   H^res = -R*T^2 * sum(z_i * d(ln phi_i)/dT)              (Eq. 9)
%   H_i^ig(T) = H_i^ig(T_ref) + integral(Cp_i dT, T_ref, T) (Eq. 10)
%   Cp_i = C1 + C2*T + C3*T^2 + C4*T^3                      (Eq. 11)
% INPUTS:
%   T          : Temperature [K]
%   P          : Pressure [Pa]
%   comp       : Mole fractions [-] (column vector, will be normalized)
%   Pc         : Critical pressures [Pa] (column vector)
%   Tc         : Critical temperatures [K] (column vector)
%   acentric   : Acentric factors [-] (column vector)
%   BIP        : Binary interaction parameter matrix [n x n]
%   Mw         : Molecular weights [g/mol] (column vector)
%   Cp_coeffs  : Ideal gas heat capacity coefficients [n x 4] matrix
%                Cp = C1 + C2*T + C3*T^2 + C4*T^3 [J/(mol·K)]
%   H_ig_ref   : Ideal gas enthalpy at T_ref=273.15 K [J/mol] (column vector)
%                (From Table 6: H_ig_ref = (H_ig/(M*R)) * M * R)
% OUTPUTS:
%   H_abs          : Absolute molar enthalpy for each component [J/mol]
%   H_mix          : Mixture absolute molar enthalpy [J/mol]
%   H_ig_specific  : Ideal gas specific enthalpy H_ig/M [J/g]
%   H_res_specific : Residual specific enthalpy H_res/M [J/g]
%   H_abs_specific : Absolute specific enthalpy H_abs/M [J/g]
% USAGE EXAMPLE (Reservoir 1 from Pedersen 2015):
%   % Table 6 values: H_ig/(M*R) in K/g
%   H_ig_per_mass_R = [-20; 20; 0; 7.5; 15; 17; 17; 25; 25; 33; ...];
%   % Convert to J/mol
%   R = 8.314462618;
%   H_ig_ref = H_ig_per_mass_R .* Mw * R;
%   
%   % Calculate enthalpies
%   [H_abs, H_mix, H_ig_spec, H_res_spec, H_abs_spec] = calculate_absolute_enthalpy(...
%       T, P, comp, Pc, Tc, acentric, BIP, Mw, Cp_coeffs, H_ig_ref);

R = 8.314462618;    % Universal gas constant [J/(mol·K)]
T_ref = 273.15;     % Reference temperature [K]


n = length(comp);
comp = comp(:);
comp = comp / sum(comp);  % Normalize mole fractions

Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Mw = Mw(:);
H_ig_ref = H_ig_ref(:);
H_ig = zeros(n, 1);

for i = 1:n
    C1 = Cp_coeffs(i, 1);
    C2 = Cp_coeffs(i, 2);
    C3 = Cp_coeffs(i, 3);
    C4 = Cp_coeffs(i, 4);
    delta_H_ig = C1 * (T - T_ref) + ...
                 C2 / 2 * (T^2 - T_ref^2) + ...
                 C3 / 3 * (T^3 - T_ref^3) + ...
                 C4 / 4 * (T^4 - T_ref^4);
    H_ig(i) = H_ig_ref(i) + delta_H_ig;
end

H_res = calculate_residual_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, R);

H_abs = H_ig + H_res;
H_mix = sum(comp .* H_abs);
H_ig_specific = H_ig ./ Mw;      % [J/g]
H_res_specific = H_res ./ Mw;    % [J/g]
H_abs_specific = H_abs ./ Mw;    % [J/g]

end

function H_res = calculate_residual_enthalpy(T, P, comp, Pc, Tc, acentric, BIP, R)
% CALCULATE_RESIDUAL_ENTHALPY - Partial molar residual enthalpies
%
% Uses central difference numerical differentiation:
%   H_res_i = -R * T^2 * (d ln(phi_i) / dT)_P,x
%           ≈ -R * T^2 * [ln(phi_i(T+dT)) - ln(phi_i(T-dT))] / (2*dT)
%
% Input/Output: All vectors are column vectors

n = length(comp);
H_res = zeros(n, 1);

% Step size for numerical differentiation
% Use relative step for better numerical stability
dT = max(0.1, 1e-4 * T);  % At least 0.1 K, or 0.01% of T

% Calculate fugacity coefficients at T-dT and T+dT
try
    [phi_minus, ~] = fugacitycoef_multicomp(comp, P, T - dT, Pc, Tc, acentric, BIP);
    [phi_plus, ~] = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP);

    % Central difference for d(ln phi)/dT
    for i = 1:n
        if phi_plus(i) > 0 && phi_minus(i) > 0
            dln_phi_dT = (log(phi_plus(i)) - log(phi_minus(i))) / (2 * dT);
            H_res(i) = -R * T^2 * dln_phi_dT;
        else
            % Fallback: use forward difference
            H_res(i) = calculate_single_component_H_res(i, T, P, comp, Pc, Tc, acentric, BIP, R, dT);
        end
    end
    
catch ME
    % If central difference fails, try forward difference for each component
    warning('calculate_absolute_enthalpy:centralDiffFailed', ...
            'Central difference failed: %s. Using forward difference.', ME.message);
    for i = 1:n
        H_res(i) = calculate_single_component_H_res(i, T, P, comp, Pc, Tc, acentric, BIP, R, dT);

    end
end

% Check for non-finite values and set to zero
bad_idx = ~isfinite(H_res);
if any(bad_idx)
    warning('calculate_absolute_enthalpy:nonFiniteValues', ...
            'Non-finite residual enthalpies detected for %d components. Setting to zero.', sum(bad_idx));
    H_res(bad_idx) = 0;
end
end

function H_res_i = calculate_single_component_H_res(i, T, P, comp, Pc, Tc, acentric, BIP, R, dT)
% Calculate residual enthalpy for a single component using forward difference

try
    [phi_0, ~] = fugacitycoef_multicomp(comp, P, T, Pc, Tc, acentric, BIP);
    [phi_forward, ~] = fugacitycoef_multicomp(comp, P, T + dT, Pc, Tc, acentric, BIP);
    
    if phi_0(i) > 0 && phi_forward(i) > 0
        dln_phi_dT = (log(phi_forward(i)) - log(phi_0(i))) / dT;
        H_res_i = -R * T^2 * dln_phi_dT;
    else
        H_res_i = 0;
    end
catch
    H_res_i = 0;
end
end