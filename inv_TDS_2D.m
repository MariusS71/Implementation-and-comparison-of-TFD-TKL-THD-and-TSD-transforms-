function [s] = inv_TDS_2D(img,S,r, c)
%% Transformata Slant 2D inversa
% INPUTS:
%  img          -- coeficientii transformarii TKL (dimensiune MxM)
%   S           -- matricea Slant ordonata  (dimensiune MxM, unde M - prima
%                 putere a lui 2 >= max(r,c) ) 
%   r           -- numarul de linii al semnalului original
%   c           -- numarul de coloane al semnalului original
%
% OUTPUT:
%   s           -- imaginea reconstruita (dimensiune rxc)
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Slant 2D inversa")

numDimensions = ndims(img);

if numDimensions == 2
  
     % Get the size of the image
    [rows, cols] = size(img);

    % Apply TDS_1D to each row
    for i = 1:rows
        img(i, :) = inv_TDS_1D(img(i, :),S);
    end

    img = img';

    for i = 1:cols
        img(i, :) = inv_TDS_1D(img(i, :),S);
    end
    
    img = img(1:r, 1:c);

    s = img;


   

elseif numDimensions == 3

    % RGB image (3D)
    
    % Extract each color channel
    img_red = img(:, :, 1);
    img_green = img(:, :, 2);
    img_blue = img(:, :, 3);

    % Get the size of the image
    [rows, cols] = size(img_red);

    % Apply TDS_1D to each row
    for i = 1:rows
        img_red(i, :) = inv_TDS_1D(img_red(i, :),S);
        img_green(i, :) = inv_TDS_1D(img_green(i, :),S);
        img_blue(i, :) = inv_TDS_1D(img_blue(i, :),S);
    end

    img_red = img_red';
    img_green = img_green';
    img_blue = img_blue';

    for i = 1:cols
        img_red(i, :) = inv_TDS_1D(img_red(i, :),S);
        img_green(i, :) = inv_TDS_1D(img_green(i, :),S);
        img_blue(i, :) = inv_TDS_1D(img_blue(i, :),S);
    end

    
    img_red = img_red(1:r, 1:c);
    img_green = img_green(1:r, 1:c);
    img_blue = img_blue(1:r, 1:c);

    s=cat(3, img_red, img_green, img_blue);
    %B=uint8(B);
    
    

else
    error('Unsupported image format. Input must be either a 2D or 3D matrix.');
end


%% SOLUTION END %%

end