%% CALCULATE DIMENSIONLESS ATTRACTION AND COVOLUME, A & B
% The Definition of Variables.
% press   : pressure
% temp    : temperature
% pressc  : critical pressure
% tempc   : critical temperature
% acentric: acentric factor
% ncomp   : the number of components

function [A, B] = calcab_multicomp(press, temp, pressc, tempc, acentric)

ncomp = size(pressc,1);
m = zeros(ncomp,1);
alpha = zeros(ncomp,1);
A = zeros(ncomp,1);
B = zeros(ncomp,1);
omegaa = 0.45724;
omegab = 0.0778;

for i = 1:ncomp
   
    pressr = press/pressc(i);
    tempr = temp/tempc(i);

    % alpha function by peng and robinson 
    % if acentric(i) > 0.49
    %     m(i) = 0.379642 + 1.48503*acentric(i) - 0.164423*acentric(i)^2 + 0.016666*acentric(i)^3;
    % else
      m(i) = 0.37464 + 1.54226*acentric(i) - 0.26992*acentric(i)^2;
   % end

    % m(i) = 0.384401 + 1.52276 * acentric(i)-0.213808* acentric(i)^2 + 0.034616 * acentric(i)^3 - 0.001976 * acentric(i)^4 ; % m as suggested by magoulus and tassios 

    alpha(i) = ( 1 + m(i)*(1 - sqrt(tempr)) )^2;
    
    
    % pina-martinez(2019) alpha rule  
    % alpha(i) = (1 + (0.3919 + 1.4996 * acentric(i) - 0.2721 * acentric(i)^2 + 0.1063*acentric(i)^3 ) * (1-tempr^(0.5))  )^2;   
    
    
    % pina-martinez(2018) alpha rule 
    % l = 0.0728 + 0.6693 * acentric(i) + 0.0925 * acentric(i)^2 ; 
    % m = 0.8788 - 0.2258* acentric(i) + 0.1695* acentric(i)^2 ;
    % alpha(i) = tempr^(2* (m-1)) * exp(l*(1-tempr^(2*m))) ; 
    

    % alpha li-yang (2012)
    % alpha(i)= exp ((0.13280-0.05052* acentric(i) + 0.25948* acentric(i)^2)* (1-tempr) + 0.81769 * log(    1+(0.31355 + 1.86745 *acentric(i)- 0.52604 * acentric(i)^2) * (1-tempr^(0.5))       )^2         ) ; 

    A(i) = omegaa*alpha(i)*pressr/tempr^2;
    B(i) = omegab*pressr/tempr;
    
end

end