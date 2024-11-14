function [s_TKL, m, Vm] = TKL_1D(s)
%% Transformata Karhunen-Loeve
% INPUTS:
%   s        -- semnnalul 1D de intrare (dimensiune N)
%
% OUTPUT:
%   s_TKL    --  vectorul coeficientilor transformatei Karhunen-Loeve
%   (Hotelling) (dimensiune N)
%   m        --  media semnalului original
%   Vm       -- vectorii proprii (dimensiune NxM)
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Karhunen-Loeve 1D")

% Ensure the signal is a column vector
s = reshape(s, [], 1);


%centrare
m = mean(s);
y = s-m;

%calcul secventa auto corelatie
disp("Calcul secventa de auto-corelatie")
r = xcov(y);
%r = xcov(y, 'unbiased');
%xcorr
r = r(length(y):end); % Extract the positive part of the auto-correlation


% matricea de auto corelatie
disp("Calcul matrice de auto-corelatie")
R=toeplitz(r);

% Calcularea valorilor și vectorilor proprii
disp("Calcul valori și vectori proprii")
[V, D] = eig(R);

% Sortarea in ordine descrescatoare
[~, order] = sort(diag(D), 'descend');
%D = D(order, order);
V = V(:, order);


% alegerea vectorilor prorii principali
Vm=V(:, 1:end); 


disp("Calcul coeficienti")
s_TKL=Vm'*y;



%% SOLUTION END %%

end