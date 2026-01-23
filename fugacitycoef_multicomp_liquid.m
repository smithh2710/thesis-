%% CALCULATE THE FUGACITY COEFFICIENT AND Z-FACTOR OF MULTI-COMPONENT SYSTEMS
% -------------------------------------------------------------------------
% The Definition of Variables.
% comp    : composition
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% BIP     : binary interaction parameter
% fugcoef : fugacity coefficient
% zfactor : compressibility factor
% -------------------------------------------------------------------------
% In this function, the minimum z-factor is automatically chosen.
function [fugcoef, zfactor] = fugacitycoef_multicomp_liquid(comp, press, temp, pressc, tempc, acentric, BIP)

[A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric);

[Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP);

zfactor = calczfactor(Amix, Bmix);

if (size(zfactor,1) > 1)
    zfactor = min(zfactor);
end

 fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2);

% fugcoef_pr = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2);
% c = [-0.2209,-0.1649,-0.1309,-0.0873,-0.0581,-0.0270]' * 1e-6;; 
% R = 8.314 ;
%     fugcoef = zeros(size(comp));
% 
%     for i = 1:length(comp)
%         fugcoef(i) = fugcoef_pr(i) / exp((c(i) * press) / (R * temp));
%     end


end
