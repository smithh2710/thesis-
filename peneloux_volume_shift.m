function [s, c, mix] = peneloux_volume_shift(Pc, Tc, acentric, comp, Vc)
% PENELOUX_VOLUME_SHIFT - Calculate Peneloux volume translation for PR-EOS
%
% Two methods available (toggle by commenting/uncommenting):
%   1. MW-based correlation (Pedersen & Christensen) - CURRENTLY ACTIVE
%   2. Rackett-based correlation - COMMENTED OUT
%
% USAGE:
%   [s, c] = peneloux_volume_shift(Pc, Tc, acentric)
%   [s, c, mix] = peneloux_volume_shift(Pc, Tc, acentric, comp, Vc)
%
% INPUTS:
%   Pc       : Critical pressures [Pa]
%   Tc       : Critical temperatures [K]
%   acentric : Acentric factors [-]
%   comp     : (Optional) Mole fractions for mixture calculation
%   Vc       : (Optional) Critical molar volumes [cm³/mol]
%
% OUTPUTS:
%   s   : Dimensionless shift parameter (s = c/b) [-]
%   c   : Volume translation parameters [m³/mol]
%   mix : (Optional) Mixture properties structure
%
% METHOD 1 - MW-BASED (Pedersen & Christensen) - CURRENTLY ACTIVE:
%   s_i = 0.0887 * ln(M_i) - 0.4668
%   c_i = s_i * b_i
%
% METHOD 2 - RACKETT-BASED - COMMENTED OUT:
%   Z_RA = 0.29056 - 0.08775 * ω
%   c_i = (0.1154 - 0.4406 * Z_RA,i) * R * Tc,i / Pc,i
%
% MIXING RULE:
%   c_mix = Σ(ζ_i * c_i)  where ζ_i = surface area fractions
%   ζ_i = (x_i * Vc_i^(2/3)) / Σ(x_j * Vc_j^(2/3))
%
% SIGN CONVENTION:
%   Positive c = volume reduction: V_corrected = V_EOS - c
%
% REFERENCES:
%   Peneloux, A., Rauzy, E., Freze, R. (1982). Fluid Phase Equilibria, 8, 7-23.
%   Pedersen, K.S., Christensen, P.L. (2007). Phase Behavior of Petroleum 
%       Reservoir Fluids, CRC Press.
%
% Author: Smit Modi
% Date: December 2025

    % Constants
    R = 8.3144598;  % [J/(mol·K)]
    Omega_b = 0.07780;  % PR-EOS covolume parameter
    
    % Molecular weights [g/mol] for Case 2 components
    % Order: C1, C2, C3, nC5, nC7, nC10
     M_gmol = [28.0 ,44.0 ,16.0 ,30.1 ,44.1 ,58.1 ,58.1 ,72.2 ,72.2 ,86.2 ,96.0 ,107 ,121 ,140.1 ,167.6, 204.7 ,243.6 ,275.3 ,317.0 ,370.4 ,456.8 ,640.9]';
    
    % Ensure column vectors
    Pc = Pc(:);
    Tc = Tc(:);
    acentric = acentric(:);
    
    nc = length(Pc);
    
    % Preallocate
    c = zeros(nc, 1);
    s = zeros(nc, 1);
    b = zeros(nc, 1);
    Z_RA = zeros(nc, 1);
    
    % Calculate for each component
    for i = 1:nc
        % Covolume b_i [m³/mol]
        b(i) = Omega_b * R * Tc(i) / Pc(i);
        
        % =====================================================================
        % METHOD 1: MW-based correlation (Pedersen & Christensen) - ACTIVE
        % =====================================================================
        s(i) = 0.0887 * log(M_gmol(i)) - 0.4668;
        c(i) = s(i) * b(i);  % [m³/mol]
        % 
        % =====================================================================
        % METHOD 2: Rackett-based correlation - COMMENTED OUT
        % Uncomment below and comment out METHOD 1 above to use Rackett
        % =====================================================================
        % Z_RA(i) = 0.29056 - 0.08775 * acentric(i);
        % c(i) = (0.1154 - 0.4406 * Z_RA(i)) * R * Tc(i) / Pc(i);
        % s(i) = c(i) / b(i);
    end
    
    % Mixture calculation (if composition provided)
    if nargin >= 4 && ~isempty(comp)
        comp = comp(:);
        comp = comp / sum(comp);  % Normalize
        
        % Calculate surface area fractions (seta)
        if nargin >= 5 && ~isempty(Vc)
            Vc = Vc(:);
            
            % Convert Vc from cm³/mol to m³/mol for internal calculation
            Vc_m3 = Vc * 1e-6;  % [m³/mol]
            
            % Surface area fractions using Vc^(2/3)
            Vc_23 = Vc_m3.^(2/3);
            sum_term = sum(comp .* Vc_23);
            
            if sum_term > 1e-15
                seta = (comp .* Vc_23) / sum_term;
            else
                seta = comp;  % Fallback to molar fractions
            end
        else
            % If no Vc provided, use molar fractions
            seta = comp;
        end
        
        % Mixture volume translation using surface area weighting
        c_mix = sum(seta .* c);  % [m³/mol]
        
        % Store mixture properties
        mix.c_mix = c_mix;           % [m³/mol]
        mix.c_mix_cm3 = c_mix * 1e6; % [cm³/mol]
        mix.seta = seta;
        mix.b = b;
        mix.s = s;
        mix.Z_RA = Z_RA;
    else
        mix = [];
    end
end