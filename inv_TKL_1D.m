function [s] = inv_TKL_1D(s_TKL, Vm, m)
%% Transformata Karhunen-Loeve inversa
% INPUTS:
%   s_TKL    --  vectorul coeficientilor transformatei Karhunen-Loeve
%   (Hotelling)
%   Vm       -- vectorii proprii 
%   m        -- media semnalului
%
% OUTPUT:
%   s        -- semnnalul 1D reconstruit
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Karhunen-Loeve 1D inversa")

% Ensure the signal is a column vector
s_TKL = reshape(s_TKL, [], 1);


% TKL invers
disp("Calcul TKL invers")
s = Vm * s_TKL + m;






%% SOLUTION END %%

end