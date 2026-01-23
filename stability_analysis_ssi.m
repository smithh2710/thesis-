function [stability, x_trial1, x_trial2, TPD_min, converged] = stability_analysis_ssi(comp, press, temp, Pc, Tc, acentric, BIP, tol, maxiter)
% PHASE STABILITY ANALYSIS USING SUCCESSIVE SUBSTITUTION ITERATION (SSI)
% Based on Michelsen (1982a,b) method from the PDF
%
% Inputs:
%   comp     : Overall composition (mole fractions)
%   press    : Pressure [Pa]
%   temp     : Temperature [K]
%   Pc       : Critical pressure vector [Pa]
%   Tc       : Critical temperature vector [K]
%   acentric : Acentric factor vector
%   BIP      : Binary interaction parameter matrix
%   tol      : Convergence tolerance (default: 1e-9)
%   maxiter  : Maximum iterations (default: 100)
%
% Outputs:
%   stability : 1 = stable (single phase), 2 = unstable (will split)
%   x_trial1  : First trial phase composition (if unstable)
%   x_trial2  : Second trial phase composition (if unstable)
%   TPD_min   : Minimum tangent plane distance value
%   converged : Convergence flag for all trials

%% Initialize parameters
if nargin < 8 || isempty(tol)
    tol = 1e-9; % From Li and Firoozabadi (2012)
end
if nargin < 9 || isempty(maxiter)
    maxiter = 100;
end

% Normalize composition
comp = comp(:) / sum(comp);
ncomp = length(comp);

% Calculate reference fugacity coefficients for feed composition z
[fugcoef_z, ~] = fugacitycoef_multicomp(comp, press, temp, Pc, Tc, acentric, BIP);
lnz_plus_lnphi_z = log(comp) + log(fugcoef_z);

%% Generate initial K-value guesses (Section 3.3.1.2 from PDF)
% For one-phase test: Nc + 4 initial guesses
K_wilson = wilsoneq(press, temp, Pc, Tc, acentric);

% Build initial K-value sets
n_initial = ncomp + 4;
K_initial_sets = zeros(ncomp, n_initial);

% First 4: Wilson-based guesses (Equation 3.43)
K_initial_sets(:, 1) = 1 ./ K_wilson;      % 1/K_wilson
K_initial_sets(:, 2) = K_wilson;           % K_wilson
K_initial_sets(:, 3) =  (K_wilson).^(1/3);   
K_initial_sets(:, 4) = (K_wilson).^(3);       

% Next Nc: Pure component guesses (Equation 3.45)
for j = 1:ncomp
    K_pure = zeros(ncomp, 1);
    K_pure(j) = 0.9 / comp(j);
    
    for i = 1:ncomp
        if i ~= j
            K_pure(i) = 0.1 / ((ncomp - 1) * comp(i));
        end
    end
    K_initial_sets(:, 4 + j) = K_pure;
end

%% Main SSI loop for all trial phases
n_trials = size(K_initial_sets, 2);
TPD_values = inf(n_trials, 1);
X_solutions = zeros(ncomp, n_trials);
converged_flags = false(n_trials, 1);

for trial = 1:n_trials
    K_initial = K_initial_sets(:, trial);
    
    % Initialize trial composition (Step 3 from PDF)
    X = comp .* K_initial;
    
    % SSI iterations (Step 4 from PDF)
    converged_trial = false;
    
    for iter = 1:maxiter
        X_old = X;
        
        % Normalize X to get mole fractions
        sum_X = sum(X);
        if sum_X > 1e-15
            x = X / sum_X;
        else
            % Near-zero composition, skip this trial
            break;
        end
        
        % Calculate fugacity coefficients for trial composition
        [fugcoef_X, ~] = fugacitycoef_multicomp(x, press, temp, Pc, Tc, acentric, BIP);
        
        % SSI update (Equation 3.42 from PDF)
        % X_i^k = exp(ln z_i + ln φ_i(z) - ln φ_i(X^(k-1)))
        for i = 1:ncomp
            X(i) = exp(lnz_plus_lnphi_z(i) - log(fugcoef_X(i)));
        end
        
        % Calculate modified TPD function (Equation 3.39)
        TPD_modified = 1.0;
        for i = 1:ncomp
            if X(i) > 1e-15
                TPD_modified = TPD_modified + X(i) * (log(X(i)) + log(fugcoef_X(i)) ...
                              - lnz_plus_lnphi_z(i) - 1);
            end
        end
        
        % Calculate beta for trivial solution check (Equation 3.49)
        beta = 0;
        for i = 1:ncomp
            beta = beta + (X(i) - comp(i)) * (log(X(i)) + log(fugcoef_X(i)) ...
                   - lnz_plus_lnphi_z(i));
        end
        
        % Calculate r parameter (Equation 3.50)
        if abs(beta) > 1e-15
            r = 2 * TPD_modified / beta;
        else
            r = 0;
        end
        
        % Check for trivial solution (Equation 3.52 - Matheis & Hickel 2017)
        if abs(r - 1) < 0.1
            % Approaching trivial solution, terminate this trial
            break;
        end
        
        % Check for divergence
        if any(X > 1e10) || any(~isfinite(X))
            % Solution diverging, terminate this trial
            break;
        end
        
        % Check convergence (Step 5.1 from PDF)
        error_norm = norm(X - X_old)^2;
        if error_norm <= tol
            converged_trial = true;
            break;
        end
    end
    
    % Store results if converged
    if converged_trial
        sum_X = sum(X);
        
        % Calculate TPD at solution point (Step 5.1)
        if sum_X > 1e-15
            TPD = -log(sum_X);           % Original TPD
            % TPD_modified = 1 - sum_X;  % Alternative (Equation 3.41)
        else
            TPD = inf;
        end
        
        X_solutions(:, trial) = X;
        TPD_values(trial) = TPD;
        converged_flags(trial) = true;
        
        % Early termination if instability detected
        if TPD < -1e-8  % Matheis and Hickel (2017) threshold
            % Unstable - no need to test other initial values
            break;
        end
    end
end

%% Determine stability from all trials
% Find minimum TPD among converged solutions
valid_indices = converged_flags & isfinite(TPD_values);

if any(valid_indices)
    valid_TPD = TPD_values(valid_indices);
    valid_X = X_solutions(:, valid_indices);
    
    [TPD_min, min_idx] = min(valid_TPD);
    
    % Normalize to get composition
    x_trial1 = valid_X(:, min_idx) / sum(valid_X(:, min_idx));
    
    % Stability criterion
    if TPD_min < -1e-8  % Instability threshold
        stability = 2;  % Unstable (will split)
        
        % Find second trial phase
        if size(valid_X, 2) > 1
            temp_TPD = valid_TPD;
            temp_TPD(min_idx) = inf;
            [~, idx_second] = min(temp_TPD);
            
            if norm(valid_X(:, idx_second) - valid_X(:, min_idx)) > 0.1
                x_trial2 = valid_X(:, idx_second) / sum(valid_X(:, idx_second));
            else
                x_trial2 = estimate_complementary_phase(x_trial1, comp);
            end
        else
            x_trial2 = estimate_complementary_phase(x_trial1, comp);
        end
    else
        stability = 1;  % Stable (single phase)
        x_trial2 = comp;
    end
    
    converged = true;
else
    % No converged solutions - assume stable
    stability = 1;
    x_trial1 = comp;
    x_trial2 = comp;
    TPD_min = 0;
    converged = false;
end

end

%% Helper function to estimate complementary phase
function x_comp = estimate_complementary_phase(x_trial, z)
% Estimate complementary phase using material balance

ncomp = length(z);
x_comp = zeros(ncomp, 1);

% Assume 50% phase split
alpha = 0.5;

for i = 1:ncomp
    % Material balance: z_i = alpha * x_trial_i + (1-alpha) * x_comp_i
    x_comp(i) = (z(i) - alpha * x_trial(i)) / (1 - alpha);
    
    % Ensure positive
    if x_comp(i) < 1e-15
        x_comp(i) = 1e-15;
    end
end

% Normalize
x_comp = x_comp / sum(x_comp);

end