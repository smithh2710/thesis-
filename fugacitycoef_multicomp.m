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
% In this function, an appropriate z-factor is automatically chosen
% according to gibbs free energy if multiple roots are found.
function [fugcoef, zfactor] = fugacitycoef_multicomp(comp, press, temp, pressc, tempc, acentric, BIP)

[A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric);

[Amix, Bmix, Amix2] = calcabmix(comp, A, B, BIP);

zfactor = calczfactor(Amix, Bmix);

if (size(zfactor,1) > 1)
    zfactor = choosezfactor(zfactor, comp, A, B, Amix, Bmix, Amix2);
end

 fugcoef = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2);

% fugcoef_pr = calcfugcoef_multicomp(zfactor, A, B, Amix, Bmix, Amix2);
% c = [-0.0718, 0.0000, -0.2209, -0.1649, -0.1309, -0.1064,-0.1064, -0.0873, -0.0873, -0.0715, -0.0217, 0.0769]' * 1e-6; 
% R = 8.314 ;
%     fugcoef = zeros(size(comp));
% 
%     for i = 1:length(comp)
%         fugcoef(i) = fugcoef_pr(i) / exp((c(i) * press) / (R * temp));
%     end
end

%% SEARTCH AND RETURN AN APPROPRIATE Z-FACTOR
% Calculate dimensionless excess gibbs free energy, and return the z
% factor which minimizes the gibbs free energy.
function minzfactor = choosezfactor(zfactor, comp, A, B, Amix, Bmix, Amix2)

gibbsenergy = [];

for i = 1:size(zfactor,1)
    fugcoef = calcfugcoef_multicomp(zfactor(i), A, B, Amix, Bmix, Amix2);
    g = calcgibbsenergy(comp, fugcoef);
    gibbsenergy = cat(1,gibbsenergy,g);
end

[~, index] = sort(gibbsenergy);
minzfactor = zfactor(index(1));

end

function g = calcgibbsenergy(comp, fugcoef)

ncomp = size(comp,1);
g = 0;

for i = 1:ncomp
    if comp(i) ~= 0
        g = g + comp(i)*log(comp(i)*fugcoef(i));
    end
end

end