function [s_slant] = TDS_1D(s,S)
%% Transformata Slant 1D
% INPUTS:
%   s          -- semnnalul 1D de intrare (dimensiune N)
%   S          -- matricea Slant ordonata  (dimensiune MxM, unde M - prima
%                 putere a lui 2 >= N)
%
% OUTPUT:
%   s_slant    -- coeficientii transformatei Slant (dimensiune M)
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Slant 1D")

% Ensure the signal is a column vector
s = reshape(s, [], 1);


%nextPowerOf2 = nextpow2(length(s));
targetLength = size(S,1);

% Zero-pad the signal if needed
if length(s) < targetLength
    s = padarray(s, [targetLength - length(s), 0], 0, 'post');
end



s_slant = S * s;
s_slant = s_slant';

%% SOLUTION END %%

end