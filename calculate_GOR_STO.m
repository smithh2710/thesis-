function [GOR, Bo, Rs] = calculate_GOR_STO(comp, P_res, T_res, Pc, Tc, acentric, BIP, M_gmol)
% CALCULATE_GOR_STO - Calculate Gas-Oil Ratio via flash to stock tank conditions
% =========================================================================
% Performs a flash calculation from reservoir conditions to stock tank
% conditions (STO: 1 atm, 15°C or 60°F) to determine:
%   - GOR: Gas-Oil Ratio [Sm³/Sm³] or [scf/STB]
%   - Bo:  Oil Formation Volume Factor [m³/Sm³]
%   - Rs:  Solution Gas-Oil Ratio [Sm³/Sm³]
%
% Method:
%   1. Flash reservoir fluid to separator conditions (optional)
%   2. Flash to stock tank conditions (1 atm, 288.15 K)
%   3. Calculate volumes of gas and oil at standard conditions
%
% INPUTS:
%   comp     : Reservoir fluid composition [mole fractions]
%   P_res    : Reservoir pressure [Pa]
%   T_res    : Reservoir temperature [K]
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   BIP      : Binary interaction parameters
%   M_gmol   : Molecular weights [g/mol]
%
% OUTPUTS:
%   GOR      : Gas-Oil Ratio [Sm³/Sm³]
%   Bo       : Oil FVF [res m³/Sm³ STO]
%   Rs       : Solution GOR [Sm³/Sm³]
%
% =========================================================================

% Standard conditions
P_STO = 101325;      % 1 atm [Pa]
T_STO = 288.15;      % 15°C [K] (or use 288.71 K for 60°F)

R = 8.3144598;       % Gas constant [J/(mol·K)]

% Ensure column vectors
comp = comp(:);
n = length(comp);

% Initialize outputs
GOR = NaN;
Bo = NaN;
Rs = NaN;

%% Step 1: Check if reservoir fluid is single phase or two-phase
% For simplicity, assume single-phase liquid in oil zone

%% Step 2: Flash to Stock Tank Conditions
try
    % Perform flash at STO conditions
    [K_sto, comp_vap_sto, comp_liq_sto, V_frac] = vaporliquideq(P_STO, T_STO, ...
        comp, Pc, Tc, acentric, BIP, 1e-10, 200);
    
    % If flash didn't converge or single phase, estimate
    if V_frac < 1e-10 || V_frac > (1-1e-10)
        % Try with Wilson K-values
        K_wilson = wilsoneq(P_STO, T_STO, Pc, Tc, acentric);
        
        % Rachford-Rice to get vapor fraction
        V_frac = solve_RR(comp, K_wilson);
        
        if V_frac < 0 || V_frac > 1
            % All liquid or all vapor
            if sum(comp .* K_wilson) > 1
                V_frac = 1;  % All vapor
            else
                V_frac = 0;  % All liquid
            end
        end
        
        comp_vap_sto = K_wilson .* comp ./ (1 + V_frac * (K_wilson - 1));
        comp_liq_sto = comp ./ (1 + V_frac * (K_wilson - 1));
    end
    
    % Normalize compositions
    comp_vap_sto = comp_vap_sto / sum(comp_vap_sto);
    comp_liq_sto = comp_liq_sto / sum(comp_liq_sto);
    
catch
    % Fallback: use Wilson K-values
    K_wilson = wilsoneq(P_STO, T_STO, Pc, Tc, acentric);
    V_frac = solve_RR(comp, K_wilson);
    
    if V_frac < 0, V_frac = 0; end
    if V_frac > 1, V_frac = 1; end
    
    comp_vap_sto = K_wilson .* comp ./ (1 + V_frac * (K_wilson - 1));
    comp_liq_sto = comp ./ (1 + V_frac * (K_wilson - 1));
    
    comp_vap_sto = comp_vap_sto / sum(comp_vap_sto);
    comp_liq_sto = comp_liq_sto / sum(comp_liq_sto);
end

%% Step 3: Calculate Molar Volumes at STO

% Gas phase - assume ideal gas at STO
V_gas_molar = R * T_STO / P_STO;  % [m³/mol]

% Liquid phase - use EOS
try
    [~, Z_liq] = fugacitycoef_multicomp_liquid(comp_liq_sto, P_STO, T_STO, ...
        Pc, Tc, acentric, BIP);
    V_liq_molar = Z_liq * R * T_STO / P_STO;  % [m³/mol]
catch
    % Fallback: estimate liquid density
    % Use average molecular weight and assume density ~700 kg/m³
    MW_liq = sum(comp_liq_sto .* M_gmol);  % [g/mol]
    rho_liq = 700;  % [kg/m³] typical for stock tank oil
    V_liq_molar = MW_liq / (rho_liq * 1000);  % [m³/mol]
end

%% Step 4: Calculate GOR

% Moles of gas and liquid per mole of feed
n_gas = V_frac;           % moles gas / mole feed
n_liq = 1 - V_frac;       % moles liquid / mole feed

% Volumes at standard conditions
V_gas_std = n_gas * V_gas_molar;  % [Sm³ gas / mol feed]
V_liq_std = n_liq * V_liq_molar;  % [Sm³ oil / mol feed]

% GOR = Volume of gas / Volume of oil at standard conditions
if V_liq_std > 1e-15
    GOR = V_gas_std / V_liq_std;  % [Sm³/Sm³]
else
    GOR = Inf;  % Dry gas
end

%% Step 5: Calculate Oil FVF (Bo)
% Bo = Volume of oil at reservoir / Volume of oil at STO

try
    [~, Z_res] = fugacitycoef_multicomp_liquid(comp, P_res, T_res, ...
        Pc, Tc, acentric, BIP);
    V_res_molar = Z_res * R * T_res / P_res;  % [m³/mol] at reservoir
catch
    % Estimate
    MW_res = sum(comp .* M_gmol);
    rho_res = 600;  % [kg/m³] typical reservoir oil
    V_res_molar = MW_res / (rho_res * 1000);
end

% Bo = V_reservoir / V_STO (per mole of liquid)
if V_liq_std > 1e-15
    Bo = V_res_molar / V_liq_std;  % [res m³ / Sm³ STO]
else
    Bo = NaN;
end

% Solution GOR (same as GOR for single-stage flash)
Rs = GOR;

end


%% Helper function: Solve Rachford-Rice equation
function V = solve_RR(z, K)
% Solve Rachford-Rice equation: sum(z_i(K_i-1)/(1+V(K_i-1))) = 0

z = z(:);
K = K(:);

% Bounds
K_max = max(K);
K_min = min(K);

if K_max < 1
    V = 0;  % All liquid
    return;
end
if K_min > 1
    V = 1;  % All vapor
    return;
end

V_min = 1 / (1 - K_max);
V_max = 1 / (1 - K_min);

V_min = max(0, V_min + 1e-10);
V_max = min(1, V_max - 1e-10);

if V_min >= V_max
    V = 0.5;
    return;
end

% Bisection
tol = 1e-12;
for iter = 1:100
    V = (V_min + V_max) / 2;
    
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