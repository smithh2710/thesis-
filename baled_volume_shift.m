function [s, c, mix] = baled_volume_shift(T, Pc, Tc, acentric, M, components, comp, Vc, verbose)
% BALED_VOLUME_SHIFT - Temperature-dependent volume translation for PR-EOS
% =========================================================================
% Reference: Baled et al. (2012), Fluid Phase Equilibria 317, 65-76
%            "Prediction of hydrocarbon densities at extreme conditions
%             using volume-translated SRK and PR equations of state fit
%             to high temperature, high pressure PVT data"
%
% Applicable range: 278-533 K, 7-276 MPa (HTHP conditions)
%
% Volume translation equation (Eq. 8):
%   c = A + B * Tr   [cm³/mol]
%
% where Tr = T/Tc is the reduced temperature
%
% For KNOWN components: A and B from Table 1 (PR column)
% For PSEUDO/UNKNOWN components: A and B from generalized correlation (Eq. 16)
%
% Generalized correlation (Eq. 16, Table 2 - PR coefficients):
%   A = k0 + k1*exp(-1/(k2*M*ω)) + k3*exp(-1/(k4*M*ω)) + k5*exp(-1/(k6*M*ω))
%   B = k0 + k1*exp(-1/(k2*M*ω)) + k3*exp(-1/(k4*M*ω)) + k5*exp(-1/(k6*M*ω))
%
% Mixture calculations:
%   seta_i = (x_i * Vc_i^(2/3)) / sum(x_i * Vc_i^(2/3))
%   Tc_mix = sum(seta_i * Tc_i)
%   A_mix = sum(x_i * A_i)  (linear mixing)
%   B_mix = sum(x_i * B_i)  (linear mixing)
%   Tr_mix = T / Tc_mix
%   c_mix = A_mix + B_mix * Tr_mix   [cm³/mol]
%
% -------------------------------------------------------------------------
% Inputs:
%   T          : Temperature [K]
%   Pc         : Critical pressures [Pa]
%   Tc         : Critical temperatures [K]
%   acentric   : Acentric factors [-]
%   M          : Molecular weights [g/mol]
%   components : Cell array of component names
%                Known: 'CO2','N2','C1','C2','C3','iC4','nC4','iC5','nC5',
%                       'nC6','nC7','nC8','iC8','nC9','nC10'
%                Use '' for pseudo-components (generalized correlation)
%   comp       : (Optional) Mole fractions [-] for mixture calculation
%   Vc         : (Optional) Critical molar volumes [cm³/mol] for mixture
%   verbose    : (Optional) Print parameter lookup table (default: false)
%
% Outputs:
%   s          : Dimensionless volume shift parameters [-] (s_i = c_i/b_i)
%   c          : Volume shift parameters [m³/mol]
%   mix        : Structure with mixture properties (if comp provided)
%                .c_mix      - mixture volume translation [m³/mol]
%                .c_mix_cm3  - mixture volume translation [cm³/mol]
%                .Tc_mix     - mixture critical temperature [K]
%                .Tr_mix     - mixture reduced temperature [-]
%                .A_mix      - mixture A parameter [cm³/mol]
%                .B_mix      - mixture B parameter [cm³/mol]
%                .A          - individual A parameters [cm³/mol]
%                .B          - individual B parameters [cm³/mol]
%                .source     - cell array indicating 'Table1' or 'Generalized'
%
% -------------------------------------------------------------------------
% Usage:
%   % Pure components
%   [s, c] = baled_volume_shift(T, Pc, Tc, acentric, M, components);
%
%   % Mixture with composition
%   [s, c, mix] = baled_volume_shift(T, Pc, Tc, acentric, M, components, comp, Vc);
%
%   % Suppress output
%   [s, c, mix] = baled_volume_shift(T, Pc, Tc, acentric, M, components, comp, Vc, false);
%
% =========================================================================

% Constants
R = 8.3144598;      % J/(mol·K)
Omega_b = 0.0778;   % PR EOS covolume constant

% Ensure column vectors
Pc = Pc(:);
Tc = Tc(:);
acentric = acentric(:);
M = M(:);

ncomp = length(Pc);

% Handle optional inputs
if nargin < 6 || isempty(components)
    components = cell(ncomp, 1);
elseif ischar(components)
    components = {components};
end

if nargin < 9 || isempty(verbose)
    verbose = false;
end

% Get Baled A and B parameters for PR-EOS
[A, B, source] = get_baled_parameters_pr(components, M, acentric, ncomp, verbose);

% Calculate reduced temperature for each component
Tr = T ./ Tc;

% Preallocate
c = zeros(ncomp, 1);
s = zeros(ncomp, 1);
c_cm3 = zeros(ncomp, 1);

for i = 1:ncomp
    % Volume shift: c = A + B * Tr [cm³/mol]
    c_cm3(i) = A(i) + B(i) * Tr(i);
    c(i) = c_cm3(i) * 1e-6;          % [m³/mol]
    
    % Covolume b_i [m³/mol]
    b_i = Omega_b * R * Tc(i) / Pc(i);
    
    % Dimensionless volume shift [-]
    s(i) = c(i) / b_i;
end

% Calculate mixture properties if composition provided
if nargin >= 7 && ~isempty(comp)
    comp = comp(:);
    comp = comp / sum(comp);  % Normalize
    
    % Linear mixing for A and B
    A_mix = sum(comp .* A);  % [cm³/mol]
    B_mix = sum(comp .* B);  % [cm³/mol]
    
    % Calculate Tc_mix using seta weighting if Vc provided
    if nargin >= 8 && ~isempty(Vc)
        Vc = Vc(:);
        Vc_23 = Vc.^(2/3);
        seta = (comp .* Vc_23) / sum(comp .* Vc_23);
        Tc_mix = sum(seta .* Tc);
        Vc_mix = sum(seta .* Vc);
        mix.seta = seta;
        mix.Vc_mix = Vc_mix;
    else
        seta = comp;
        Tc_mix = sum(comp .* Tc);
        mix.seta = seta;
        mix.Vc_mix = [];
    end
    
    % Mixture reduced temperature
    Tr_mix = T / Tc_mix;
    
    % Mixture volume translation
    c_mix_cm3 = A_mix + B_mix * Tr_mix;  % [cm³/mol]
    c_mix = c_mix_cm3 * 1e-6;            % [m³/mol]
    
    % Store mixture properties
    mix.Tc_mix = Tc_mix;
    mix.Tr_mix = Tr_mix;
    mix.A_mix = A_mix;
    mix.B_mix = B_mix;
    mix.c_mix = c_mix;
    mix.c_mix_cm3 = c_mix_cm3;
    mix.A = A;
    mix.B = B;
    mix.c_i_cm3 = c_cm3;
    mix.source = source;
else
    mix = [];
end

end


function [A, B, source] = get_baled_parameters_pr(components, M, acentric, ncomp, verbose)
% GET_BALED_PARAMETERS_PR - Get Baled A and B parameters for PR EOS
%
% Table 1 values for known components (PR column)
% Eq. 16 with Table 2 coefficients for pseudo-components
%
% =========================================================================
% Table 1: PR-EOS Volume Translation Parameters from Baled et al. (2012)
% c = A + B * Tr  [cm³/mol]
% =========================================================================

% Database: {Name, A [cm³/mol], B [cm³/mol]}
table1 = {
    % Light gases
    'N2',         -1.638,    -0.698      % Nitrogen
    'nitrogen',   -1.638,    -0.698
    'CO2',        -1.181,     0.756      % Carbon dioxide
    'carbon dioxide', -1.181, 0.756
    
    % n-Alkanes
    'C1',         -3.047,    -0.610      % Methane
    'methane',    -3.047,    -0.610
    'C2',          0.364,    -0.513      % Ethane
    'ethane',      0.364,    -0.513
    'C3',         -3.328,    -3.189      % Propane
    'propane',    -3.328,    -3.189
    'nC4',         2.590,   -10.64       % n-Butane
    'n-butane',    2.590,   -10.64
    'C4',          2.590,   -10.64
    'nC5',         7.181,   -13.89       % n-Pentane
    'n-pentane',   7.181,   -13.89
    'C5',          7.181,   -13.89
    'nC6',         8.940,   -14.32       % n-Hexane
    'n-hexane',    8.940,   -14.32
    'C6',          8.940,   -14.32
    'nC7',        11.24,    -14.57       % n-Heptane
    'n-heptane',  11.24,    -14.57
    'C7',         11.24,    -14.57
    'nC8',        20.70,    -23.73       % n-Octane
    'n-octane',   20.70,    -23.73
    'C8',         20.70,    -23.73
    'nC9',        25.59,    -25.38       % n-Nonane
    'n-nonane',   25.59,    -25.38
    'C9',         25.59,    -25.38
    'nC10',       33.71,    -30.91       % n-Decane
    'n-decane',   33.71,    -30.91
    'C10',        33.71,    -30.91
    
    % Branched alkanes
    'iC4',         1.490,    -9.540      % Isobutane
    'i-butane',    1.490,    -9.540
    'isobutane',   1.490,    -9.540
    'iC5',         5.752,   -12.60       % Isopentane
    'i-pentane',   5.752,   -12.60
    'isopentane',  5.752,   -12.60
    'iC8',         7.824,   -19.51       % Isooctane (2,2,4-trimethylpentane)
    'isooctane',   7.824,   -19.51
};

% =========================================================================
% Table 2: Generalized Correlation Coefficients (Eq. 16) for PR-EOS
% A = k0 + k1*exp(-1/(k2*M*ω)) + k3*exp(-1/(k4*M*ω)) + k5*exp(-1/(k6*M*ω))
% B = k0 + k1*exp(-1/(k2*M*ω)) + k3*exp(-1/(k4*M*ω)) + k5*exp(-1/(k6*M*ω))
% =========================================================================
kA = [-4.1034, 31.723, 0.0531, 188.68, 0.0057, 20196, 0.0003];
kB = [-0.3489, -28.547, 0.0687, -817.73, 0.0007, -65.067, 0.0076];

% Initialize outputs
A = zeros(ncomp, 1);
B = zeros(ncomp, 1);
source = cell(ncomp, 1);

if verbose
    fprintf('\n--- Baled Volume Shift: Parameter Lookup ---\n');
    fprintf('%-12s %12s %12s %12s\n', 'Component', 'A [cm³/mol]', 'B [cm³/mol]', 'Source');
    fprintf('%s\n', repmat('-', 1, 52));
end

for i = 1:ncomp
    found = false;
    comp_name = '';
    
    % Check if component name is provided
    if i <= length(components) && ~isempty(components{i})
        comp_name = strtrim(components{i});  % Remove leading/trailing spaces
        
        % Search in Table 1
        for j = 1:size(table1, 1)
            if strcmpi(comp_name, table1{j, 1})
                A(i) = table1{j, 2};
                B(i) = table1{j, 3};
                source{i} = 'Table1';
                found = true;
                break;
            end
        end
    end
    
    % Use generalized correlation if not found (Eq. 16)
    if ~found
        Mw = M(i) * acentric(i);
        
        % Prevent division by zero or very small values
        if Mw < 0.01
            Mw = 0.01;
        end
        
        % Calculate A using Eq. 16
        A(i) = kA(1) + kA(2)*exp(-1/(kA(3)*Mw)) + ...
               kA(4)*exp(-1/(kA(5)*Mw)) + kA(6)*exp(-1/(kA(7)*Mw));
        
        % Calculate B using Eq. 16
        B(i) = kB(1) + kB(2)*exp(-1/(kB(3)*Mw)) + ...
               kB(4)*exp(-1/(kB(5)*Mw)) + kB(6)*exp(-1/(kB(7)*Mw));
        
        source{i} = 'Generalized';
        
        if isempty(comp_name)
            comp_name = sprintf('Pseudo%d', i);
        end
    end
    
    if verbose
        fprintf('%-12s %12.4f %12.4f %12s\n', comp_name, A(i), B(i), source{i});
    end
end

if verbose
    fprintf('%s\n', repmat('-', 1, 52));
end

end