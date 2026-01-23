function [s, c, mix] = abudour_volume_shift(comp, press, temp, Pc, Tc, acentric, Vc, components, BIP, verbose)
% ABUDOUR_VOLUME_SHIFT - Volume translation for PR-EOS using Abudour et al. (2012, 2013)
% =========================================================================
% Reference: 
%   Abudour et al. (2012), Fluid Phase Equilibria 335, 74-87 (Pure fluids)
%   Abudour et al. (2013), Fluid Phase Equilibria 349, 37-55 (Mixtures)
%
% This function uses:
%   - abudour_parameters.m for KNOWN components (optimized c1 from Table 1)
%   - abudour_generalized_c1.m for PSEUDO/UNKNOWN components (Eq. 11)
%
% Key Equations:
%   V_VTPR = V_PR + c - δc * (0.35 / (0.35 + d))           (Eq. 7)
%   δc = (R*Tc/Pc) * (Zc_EOS - Zc_exp)                     (Eq. 8)
%   c = (R*Tc/Pc) * [c1 - (0.004 + c1) * exp(-2*d)]        (Eq. 9)
%   d = (1/(R*Tc)) * |∂P/∂ρ|_T                             (Eq. 6)
%   c1 = 0.4266 * Zc - 0.1101  (generalized, Eq. 11)
%
% -------------------------------------------------------------------------
% Inputs:
%   comp            : Mole fractions [-] (column vector)
%   press           : Pressure [Pa]
%   temp            : Temperature [K]
%   Pc              : Critical pressures [Pa] (column vector)
%   Tc              : Critical temperatures [K] (column vector)
%   acentric        : Acentric factors [-] (column vector)
%   Vc              : Critical molar volumes [cm³/mol] (column vector)
%   components      : Cell array of component names (e.g., {'CO2','C1','C2',''})
%                     Use '' for pseudo-components (will use generalized c1)
%   BIP             : Binary interaction parameter matrix (optional)
%   verbose         : Print parameter lookup table (default: false)
%
% Outputs:
%   s   : Dimensionless volume shift parameters [-] for each component
%   c   : Volume shift parameters [m³/mol] for each component
%   mix : Structure with mixture properties:
%         .c_mix      - effective mixture volume translation [m³/mol]
%         .c_mix_cm3  - effective mixture volume translation [cm³/mol]
%         .delta_c    - critical point correction [m³/mol]
%         .d          - distance function [-]
%         .c1_mix     - mixture c1 parameter [-]
%         .c1         - individual c1 parameters [-]
%         .Tc_mix     - mixture critical temperature [K]
%         .Pc_mix     - mixture critical pressure [Pa]
%         .Vc_mix     - mixture critical volume [cm³/mol]
%         .V_PR       - untranslated PR volume [m³/mol]
%         .V_VTPR     - translated volume [m³/mol]
%         .source     - cell array indicating 'Table1' or 'Generalized'
%
% -------------------------------------------------------------------------
% Usage:
%   % With component names (uses optimized c1 when available)
%   names = {'CO2', 'N2', 'C1', 'C2', 'C3', '', ''};  % '' for pseudo-components
%   [s, c, mix] = abudour_volume_shift(comp, P, T, Pc, Tc, w, Vc, names, BIP);
%
%   % Without names (uses generalized correlation for all)
%   [s, c, mix] = abudour_volume_shift(comp, P, T, Pc, Tc, w, Vc, {}, BIP);
%
%   % With verbose output
%   [s, c, mix] = abudour_volume_shift(comp, P, T, Pc, Tc, w, Vc, names, BIP, true);
%
% =========================================================================

R = 8.3144598;           % J/(mol·K)
Zc_EOS = 0.3074;         % PR-EOS critical compressibility
Omega_b = 0.0778;        % PR-EOS b parameter constant

% Ensure column vectors
nc = length(comp);
comp = comp(:);
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
Vc = Vc(:);

% Handle optional component_names input
if nargin < 8 || isempty(components)
    components = cell(nc, 1);
end

% Handle optional BIP input
if nargin < 9 || isempty(BIP)
    BIP = zeros(nc, nc);
end

% Handle verbose flag
if nargin < 10 || isempty(verbose)
    verbose = false;
end

%% Get c1 and Zc for each component
c1 = zeros(nc, 1);
Zc = zeros(nc, 1);
source = cell(nc, 1);

% if verbose
%     fprintf('\n--- Abudour Volume Shift: c1 Parameter Lookup ---\n');
%     fprintf('%-12s %12s %10s %12s\n', 'Component', 'c1', 'Zc', 'Source');
%     fprintf('%s\n', repmat('-', 1, 50));
% end

for i = 1:nc
    comp_name = '';
    
    % Check if component name is provided
    if i <= length(components) && ~isempty(components{i})
        comp_name = strtrim(components{i});
        
        % Try to get from database
        [c1_db, zc_db, ~, ~, ~, ~, found] = abudour_parameters(comp_name);
        
        if found
            c1(i) = c1_db;
            Zc(i) = zc_db;
            source{i} = 'Table1';
        else
            % Use generalized correlation
            [c1(i), Zc(i)] = abudour_generalized_c1(acentric(i));
            source{i} = 'Generalized';
        end
    else
        % Pseudo-component: use generalized correlation
        [c1(i), Zc(i)] = abudour_generalized_c1(acentric(i));
        source{i} = 'Generalized';
        comp_name = sprintf('Pseudo%d', i);
    end
    
    if verbose
        fprintf('%-12s %12.6f %10.4f %12s\n', comp_name, c1(i), Zc(i), source{i});
    end
end

if verbose
    fprintf('%s\n', repmat('-', 1, 50));
end

%% Calculate EOS parameters for each component
a_coef = zeros(nc, 1);
b_coef = zeros(nc, 1);
for i = 1:nc
    [a_coef(i), b_coef(i)] = calc_ab_PR(Tc(i), Pc(i), acentric(i), temp);
end

%% Calculate mixture EOS parameters using classical mixing rules
a_mix = 0;
b_mix = 0;
for i = 1:nc
    for j = 1:nc
        a_ij = sqrt(a_coef(i) * a_coef(j)) * (1 - BIP(i,j));
        a_mix = a_mix + comp(i) * comp(j) * a_ij;
        b_ij = (b_coef(i) + b_coef(j)) / 2;
        b_mix = b_mix + comp(i) * comp(j) * b_ij;
    end
end

%% Calculate mixture critical properties (Chueh-Prausnitz correlations)
% Surface fraction (Eq. 22)
seta = (comp .* Vc.^(2/3)) / sum(comp .* Vc.^(2/3));

% Mixture critical temperature (Eq. 23)
Tc_mix = sum(seta .* Tc);

% Mixture critical volume (Eq. 21)
Vc_mix = sum(seta .* Vc);  % [cm³/mol]

% Mixture acentric factor (Eq. 25)
omega_mix = sum(comp .* acentric);

% Mixture critical pressure (Eq. 24 - Aalto et al.)
Vc_mix_m3 = Vc_mix * 1e-6;  % [m³/mol]
Pc_mix = (0.2905 - 0.085 * omega_mix) * R * Tc_mix / Vc_mix_m3;  % [Pa]

% Mixture Zc
Zc_mix = sum(comp .* Zc);

%% Calculate mixture c1 (Eq. 17 - linear mixing rule)
c1_mix = sum(comp .* c1);

%% Calculate molar volume from untranslated PR-EOS
[V_PR, ~] = calc_PR_volume(press, temp, a_mix, b_mix, R);

%% Calculate distance function (Eq. 6, 18 for mixtures)
d = calc_distance_function(V_PR, temp, a_mix, b_mix, R, Tc_mix);

%% Calculate volume translation for mixture (Eq. 9, 16)
c_mix = (R * Tc_mix / Pc_mix) * (c1_mix - (0.004 + c1_mix) * exp(-2 * d));  % [m³/mol]

%% Calculate delta_c for mixture (Eq. 8, 19)
V_PR_cm = (R * Tc_mix / Pc_mix) * Zc_EOS;  % [m³/mol] - PR predicted critical volume
V_cm = Vc_mix * 1e-6;                       % [m³/mol] - true critical volume
delta_c = V_PR_cm - V_cm;                   % [m³/mol]

%% Final effective volume translation (Eq. 7, 15)
c_eff = - ( c_mix - delta_c * (0.35 / (0.35 + d)));  % [m³/mol]

%% Calculate individual component volume shifts
c = zeros(nc, 1);
s = zeros(nc, 1);
for i = 1:nc
    % Individual component c using their own properties
    [V_i, ~] = calc_PR_volume(press, temp, a_coef(i), b_coef(i), R);
    d_i = calc_distance_function(V_i, temp, a_coef(i), b_coef(i), R, Tc(i));
    
    % Volume translation (Eq. 9)
    c_i = (R * Tc(i) / Pc(i)) * (c1(i) - (0.004 + c1(i)) * exp(-2 * d_i));
    
    % Critical point correction (Eq. 8)
    delta_c_i = (R * Tc(i) / Pc(i)) * (Zc_EOS - Zc(i));
    
    % Effective volume shift (Eq. 7)
    c(i) = - ( c_i - delta_c_i * (0.35 / (0.35 + d_i)) ) ;  % [m³/mol]
    
    % Dimensionless shift
    s(i) = c(i) / b_coef(i);
end

%% Output structure
mix.c_mix = c_eff;                    % Effective mixture volume translation [m³/mol]
mix.c_mix_cm3 = c_eff * 1e6;          % [cm³/mol]
mix.c_before_correction = c_mix;      % c before δc correction [m³/mol]
mix.delta_c = delta_c;                % Critical point correction [m³/mol]
mix.d = d;                            % Distance function [-]
mix.c1_mix = c1_mix;                  % Mixture c1 parameter [-]
mix.c1 = c1;                          % Individual c1 parameters [-]
mix.Zc = Zc;                          % Individual Zc values [-]
mix.Tc_mix = Tc_mix;                  % Mixture critical temperature [K]
mix.Pc_mix = Pc_mix;                  % Mixture critical pressure [Pa]
mix.Vc_mix = Vc_mix;                  % Mixture critical volume [cm³/mol]
mix.Zc_mix = Zc_mix;                  % Mixture critical compressibility [-]
mix.omega_mix = omega_mix;            % Mixture acentric factor [-]
mix.V_PR = V_PR;                      % Untranslated PR volume [m³/mol]
mix.V_VTPR = V_PR - c_eff;            % Translated volume [m³/mol]
mix.seta = seta;                      % Surface fractions [-]
mix.source = source;                  % Source of c1 ('Table1' or 'Generalized')

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [a, b] = calc_ab_PR(Tc, Pc, omega, T)
% Calculate PR-EOS a and b parameters with alpha function
% Using standard PR alpha function

R = 8.3144598;

% PR-EOS constants
Omega_a = 0.45724;
Omega_b = 0.07780;

% Calculate b (temperature independent)
b = Omega_b * R * Tc / Pc;

% Calculate m for alpha function
if omega > 0.49
    m = 0.379642 + 1.48503*omega - 0.164423*omega^2 + 0.016666*omega^3;
else
    m = 0.37464 + 1.54226*omega - 0.26992*omega^2;
end

% Calculate alpha
Tr = T / Tc;
alpha = (1 + m * (1 - sqrt(Tr)))^2;

% Calculate a
a = Omega_a * (R * Tc)^2 / Pc * alpha;

end

function [V, rho] = calc_PR_volume(P, T, a, b, R)
% Solve PR-EOS for molar volume (liquid root)

% Convert to cubic in Z
A = a * P / (R * T)^2;
B = b * P / (R * T);

p2 = -(1 - B);
p1 = A - 3*B^2 - 2*B;
p0 = -(A*B - B^2 - B^3);

% Solve cubic
roots_Z = roots([1, p2, p1, p0]);

% Get real roots
real_roots = roots_Z(abs(imag(roots_Z)) < 1e-10);
real_roots = real(real_roots);

% For liquid, take smallest positive root > B
valid_roots = real_roots(real_roots > B);
if isempty(valid_roots)
    Z = max(real_roots);
else
    Z = min(valid_roots);  % Liquid root
end

V = Z * R * T / P;  % [m³/mol]
rho = 1 / V;        % [mol/m³]

end

function d = calc_distance_function(V, T, a, b, R, Tc)
% Calculate dimensionless distance function (Eq. 6)
% d = (1/(R*Tc)) * |∂P/∂ρ|_T

% ∂P/∂V for PR-EOS
denom1 = V - b;
denom2 = V*(V + b) + b*(V - b);  % = V² + 2bV - b²

dPdV = -R*T / denom1^2 + a * (2*V + 2*b) / denom2^2;

% ∂P/∂ρ = -V² * ∂P/∂V
dPdrho = -V^2 * dPdV;

% Distance function
d = abs(dPdrho) / (R * Tc);

end