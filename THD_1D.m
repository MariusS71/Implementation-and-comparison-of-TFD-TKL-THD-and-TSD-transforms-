function [s_THD, N] = THD_1D(s)
%% Transformata Hartley 1D
% INPUTS:
%   s        -- semnnalul 1D de intrare   
%
% OUTPUT:
%   s_THD    --  vectorul coeficientilor transformatei Hartley
%   N        --  dimensiunea semnalului original
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Hartley 1D")

% Ensure the signal is a column vector
s = reshape(s, [], 1);

N=size(s,1);


%formez matricea de transformare pentru forma matriceala
disp("Calcul forma matriceala")
n = 0:(N-1);
k=n;
[n,k] = meshgrid(n,k);
H= cos(2*pi*n.*k/N) + sin(2*pi*n.*k/N);

%obtinerea coeficientilor
disp("Obtinere coeficienti")
s_THD = H*s;




%% SOLUTION END %%

end