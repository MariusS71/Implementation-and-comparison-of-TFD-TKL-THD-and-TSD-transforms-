function y = inv_THD_1D(s_THD, N)
%% Transformata Hartley 1D inversa
% INPUTS:
%   s_THD       --  vectorul coeficientilor transformatei Hartley
%   N           -- dimensiunea semnalului original
%
% OUTPUT:
%   y    --  semnnalul 1D reconstruit
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Hartley inversa 1D")

% Ensure the signal is a column vector
s_THD = reshape(s_THD, [], 1);


%formez matricea de transformare pentru forma matriceala
disp("Calcul forma matriceala")
n = 0:(N-1);
k=n;
[n,k] = meshgrid(n,k);
H= cos(2*pi*n.*k/N) + sin(2*pi*n.*k/N);


% aplicare inversa pentru obtinerea semnalului original
disp("Obtinere semnal original")
y =  (H * s_THD) / N;




%% SOLUTION END %%

end