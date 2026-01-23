

%% Input Data (from PDF Table 3.2)

components = {'C3H8', 'n-C4H10'};
Tc = [369.83; 425.12];      % Critical temperature [K]  
Pc = [42.48; 37.96] * 1e5;  % Critical pressure [Pa]
omega = [0.1523; 0.2002];   % Acentric factors

% Binary interaction parameters
BIP = [0.0000, 0.0026;      
       0.0026, 0.0000];

% System conditions
temp = 343.17;               % Temperature [K]
press = 14 * 1e5;           % Pressure [Pa] (14 bar)
z = [0.5; 0.5];             % Feed composition

% Tolerance and max iterations
tol = 1e-9;
maxiter = 100;


[stability, x_trial1, x_trial2, TPD_min, ~] = stability_analysis_ssi(z, press, temp, Pc, Tc, omega, BIP, tol, maxiter)

%[stability, x_trial1, x_trial2, TPD_min] = stability_analysis_qnss(z, press, temp, Pc, Tc, omega, BIP, tol, maxiter)