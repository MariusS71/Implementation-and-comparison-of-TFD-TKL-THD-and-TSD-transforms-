function [s] = inv_TDS_1D(s_slant,S)
%% Transformata Slant 1D inversa
% INPUTS:
%   s_slant     -- coeficientii transformatei Slant (dimensiune M)
%   S           -- matricea Slant ordonata  (dimensiune MxM)
%
% OUTPUT:
%   s           --  semnalul original reconstruit (dimensiune M)
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Slant 1D inversa")

% Ensure the signal is a column vector
s_slant = reshape(s_slant, [], 1);


s =  s_slant' * S;
s=s';



%% SOLUTION END %%

end