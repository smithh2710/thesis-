% ENHANCED STABILITY ANALYSIS USING TANGENT PLANE DISTANCE (TPD)
% Based on Michelsen (1982) and Nghiem & Li (1984) QNSS method

% Inputs:
%   comp     : Overall composition (mole fractions)
%   press    : Pressure [Pa]
%   temp     : Temperature [K]
%   Pc       : Critical pressure vector [Pa]
%   Tc       : Critical temperature vector [K]
%   acentric : Acentric factor vector
%   BIP      : Binary interaction parameter matrix
%   tol      : Convergence tolerance (default: 1e-8)
%   maxiter  : Maximum iterations (default: 100)
%
% Outputs:
%   stability : 1 = stable (single phase), 2 = unstable (will split)
%   x_trial1  : First trial phase composition (if unstable)
%   x_trial2  : Second trial phase composition (if unstable)
%   TPD_min   : Minimum tangent plane distance value



function [stability, x_trial1, x_trial2, TPD_min] = stability_analysis_qnss(comp, press, temp, Pc, Tc, acentric, BIP, tol, maxiter)


% Set default values if not provided
if nargin < 8 || isempty(tol)
    tol = 1e-8;
end
if nargin < 9 || isempty(maxiter)
    maxiter = 100;
end

% Normalize composition
comp = comp(:) / sum(comp);
ncomp = length(comp);

% Calculate reference fugacity coefficients for overall composition
[fugcoef_ref, ~] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP);
lnfugcoef_ref = log(fugcoef_ref);

% Initialize storage for multiple trial calculations
n_trials = 5;  % 4 Wilson-based + 1 pure component trial
TPD = zeros(n_trials, 1);
u_min = zeros(n_trials, ncomp);
x_min = zeros(n_trials, ncomp);

%% TRIAL 1-4: Wilson equation based initializations
for trial = 1:4
    % Different initialization strategies
    switch trial
        case 1  % Vapor-like initialization
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, 1.0);
        case 2  % Intermediate vapor
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, 0.33);
        case 3  % Liquid-like initialization
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, -1.0);
        case 4  % Intermediate liquid
            lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, -0.33);
    end
    
    % Perform stability test for both vapor-like and liquid-like phases
    [u_vapor, converged_v] = qnss_iteration(comp, lnK, press, temp, Pc, Tc,acentric, BIP, lnfugcoef_ref, tol, maxiter, 'vapor');
    [u_liquid, converged_l] = qnss_iteration(comp, lnK, press, temp, Pc, Tc,acentric, BIP, lnfugcoef_ref, tol, maxiter, 'liquid');
    
    % Select the trial with larger sum(u)
    if sum(u_vapor) > sum(u_liquid)
        u_min(trial, :) = u_vapor;
    else
        u_min(trial, :) = u_liquid;
    end
    
    % Calculate TPD for this trial
    x_min(trial, :) = u_min(trial, :) / sum(u_min(trial, :));
    [fugcoef_trial, ~] = fugacitycoef_multicomp(x_min(trial, :)', press, temp, Pc, Tc, acentric, BIP);
    
    TPD(trial) = calculate_tpd(x_min(trial, :), comp', fugcoef_trial, fugcoef_ref);
end

%% TRIAL 5: Pure component initialization (most volatile component)
% Find the most volatile component (lowest boiling point approximation)
[~, most_volatile] = min(Tc .* sqrt(Pc));

lnK = zeros(ncomp, 1);
for i = 1:ncomp
    if i == most_volatile
        lnK(i) = log((1 - 1e-15) / comp(i));
    else
        lnK(i) = log(1e-15 / (ncomp - 1) / comp(i));
    end
end

% Perform stability test
[u_vapor, ~] = qnss_iteration(comp, lnK, press, temp, Pc, Tc, ...
                             acentric, BIP, lnfugcoef_ref, tol, maxiter, 'vapor');
[u_liquid, ~] = qnss_iteration(comp, lnK, press, temp, Pc, Tc, ...
                              acentric, BIP, lnfugcoef_ref, tol, maxiter, 'liquid');

if sum(u_vapor) > sum(u_liquid)
    u_min(5, :) = u_vapor;
else
    u_min(5, :) = u_liquid;
end

x_min(5, :) = u_min(5, :) / sum(u_min(5, :));
[fugcoef_trial, ~] = fugacitycoef_multicomp(x_min(5, :)', press, temp, Pc, Tc, acentric, BIP);
TPD(5) = calculate_tpd(x_min(5, :), comp', fugcoef_trial, fugcoef_ref);

%% Determine stability and output results
[TPD_min, idx_min] = min(TPD);

% Find a second distinct minimum if it exists
TPD_temp = TPD;
TPD_temp(idx_min) = inf;
[TPD_second, idx_second] = min(TPD_temp);

% Stability criterion: TPD < -tolerance indicates instability
if TPD_min < -1e-6
    stability = 2;  % Unstable (will split into phases)
    x_trial1 = x_min(idx_min, :)';
    
    % Find second trial phase that is sufficiently different
    if abs(TPD_second - TPD_min) < 0.001 || norm(x_min(idx_second, :) - x_min(idx_min, :)) < 0.1
        % Use complementary phase estimate
        x_trial2 = estimate_complementary_phase(x_trial1, comp);
    else
        x_trial2 = x_min(idx_second, :)';
    end
else
    stability = 1;  % Stable (single phase)
    x_trial1 = comp;
    x_trial2 = comp;
end

end

%% QNSS ITERATION FUNCTION
function [u, converged] = qnss_iteration(z, lnK_init, press, temp, Pc, Tc, acentric, BIP, lnfugcoef_z, tol, maxiter, phase_type)
% Quasi-Newton Successive Substitution iteration for stability analysis

ncomp = length(z);
lnK = lnK_init(:);

% Initial guess based on phase type
if strcmp(phase_type, 'vapor')
    u = exp(lnK) .* z;
else
    u = z ./ exp(lnK);
end

% Normalize
x = u / sum(u);

% Calculate initial fugacity coefficients
[fugcoef_x, ~] = fugacitycoef_multicomp(x, press, temp, Pc, Tc, acentric, BIP);
lnfugcoef_x = log(fugcoef_x);

% Initial residual
d = log(u) + lnfugcoef_x - log(z) - lnfugcoef_z;
err = norm(d);

% QNSS parameters
sigma = 1.0;
iter = 0;
converged = false;

while err > tol && iter < maxiter
    % Store old values
    d_old = d;
    
    % Update step
    dlnK = -sigma * d;
    
    % Limit step size
    max_step = 3.0;
    dlnK = sign(dlnK) .* min(abs(dlnK), max_step);
    
    % Update lnK
    lnK = lnK + dlnK;
    
    % Update composition
    u = exp(lnK) .* z;
    x = u / sum(u);
    
    % Calculate new fugacity coefficients
    [fugcoef_x, ~] = fugacitycoef_multicomp(x, press, temp, Pc, Tc, acentric, BIP);
    lnfugcoef_x = log(fugcoef_x);
    
    % New residual
    d = log(u) + lnfugcoef_x - log(z) - lnfugcoef_z;
    
    % Update sigma using quasi-Newton formula
    dd = d - d_old;
    denominator = dlnK' * dd;
    if abs(denominator) > 1e-10
        sigma_new = abs(-(dlnK' * d_old) / denominator) * sigma;
        sigma = min(max(sigma_new, 0.1), 2.0);  % Limit sigma range
    end
    
    err = norm(d);
    iter = iter + 1;
    
    % Reset sigma periodically
    if mod(iter, 10) == 0
        sigma = 1.0;
    end
end

if iter < maxiter
    converged = true;
end

end

%% WILSON K-FACTOR CALCULATION
function lnK = wilson_K_factors(press, temp, Pc, Tc, acentric, factor)
% Calculate Wilson equation K-factors with scaling factor

ncomp = length(Pc);
lnK = zeros(ncomp, 1);

for i = 1:ncomp
    lnK(i) = factor * (5.373 * (1 + acentric(i)) * (1 - Tc(i)/temp) + log(Pc(i)/press));
end

end

%% TANGENT PLANE DISTANCE CALCULATION
function TPD = calculate_tpd(x_trial, z, fugcoef_trial, fugcoef_z)
% Calculate tangent plane distance

ncomp = length(z);
TPD = 0;

for i = 1:ncomp
    if x_trial(i) > 1e-15 && z(i) > 1e-15
        TPD = TPD + x_trial(i) * (log(x_trial(i)) + log(fugcoef_trial(i)) - ...
                                  log(z(i)) - log(fugcoef_z(i)));
    end
end

end

%% ESTIMATE COMPLEMENTARY PHASE
function x_comp = estimate_complementary_phase(x_trial, z)
% Estimate a complementary phase composition using lever rule approximation

ncomp = length(z);
x_comp = zeros(ncomp, 1);

% Simple complementary estimate
alpha = 0.5;  % Phase fraction guess
for i = 1:ncomp
    x_comp(i) = (z(i) - alpha * x_trial(i)) / (1 - alpha);
    x_comp(i) = max(x_comp(i), 1e-15);  % Ensure positive
end

% Normalize
x_comp = x_comp / sum(x_comp);

end