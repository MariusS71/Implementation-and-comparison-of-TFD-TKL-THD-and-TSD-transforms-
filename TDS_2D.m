function [s_slant] = TDS_2D(img, S)
%% Transformata Slant 2D
% INPUTS:
%   img         -- semnnalul 2D de intrare (dimensiune N1xN2)
%   S           -- matricea Slant ordonata  (dimensiune MxM, unde M - prima
%                 putere a lui 2 >= max(N1,N2) )
%
% OUTPUT:
%   s_TDS       -- coeficientii transformatei Slant (MxM)
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Slant 2D")

numDimensions = ndims(img);

if numDimensions == 2
   

    % Get the size of the image
    [rows, cols] = size(img);


    % Find the next power of two
    nextPowerOf2 = nextpow2(max(cols,rows));
    targetLength = 2^nextPowerOf2;

    %img = padarray(img, [targetLength - rows, targetLength - cols], 0, 'post');


    aux = zeros(targetLength,rows);
    aux2= zeros(targetLength, targetLength);

    % Apply TDS_1D to each row
    for i = 1:rows
        aux(:, i) = TDS_1D(img(i, :), S)';
    end

   % img = img';

    for i = 1:targetLength
        aux2(i, :) = TDS_1D(aux(i, :), S);
    end



    s_slant = aux2;


elseif numDimensions == 3
    % RGB image (3D)
    
    % Extract each color channel
    img_red = img(:, :, 1);
    img_green = img(:, :, 2);
    img_blue = img(:, :, 3);


    
    % Get the size of the image
    [rows, cols] = size(img_red);

    % Find the next power of two
    nextPowerOf2 = nextpow2(max(cols,rows));
    targetLength = 2^nextPowerOf2;

    %img_red = padarray(img_red, [targetLength - rows, targetLength - cols], 0, 'post');
    %img_green = padarray(img_green, [targetLength - rows, targetLength - cols], 0, 'post');
    %img_blue = padarray(img_blue, [targetLength - rows, targetLength - cols], 0, 'post');


    % Apply TDS_1D to each row
    aux_red = zeros(targetLength,rows);
    aux_green = zeros(targetLength,rows);
    aux_blue = zeros(targetLength,rows);

    aux2_red= zeros(targetLength, targetLength);
    aux2_green= zeros(targetLength, targetLength);
    aux2_blue= zeros(targetLength, targetLength);


    % Apply TDS_1D to each row
    for i = 1:rows
        aux_red(:, i) = TDS_1D(img_red(i, :), S)';
        aux_green(:, i) = TDS_1D(img_green(i, :), S)';
        aux_blue(:, i) = TDS_1D(img_blue(i, :), S)';
    end


    for i = 1:targetLength
        aux2_red(i, :) = TDS_1D(aux_red(i, :), S);
        aux2_green(i, :) = TDS_1D(aux_green(i, :), S);
        aux2_blue(i, :) = TDS_1D(aux_blue(i, :), S);
    end


    s_slant =  cat(3, aux2_red, aux2_green, aux2_blue);
   

else
    error('Unsupported image format. Input must be either a 2D or 3D matrix.');
end


%% SOLUTION END %%

end