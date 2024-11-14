function [img_THD, NL, NC] = THD_2D(img)
%% Transformata Hartley 2D
% INPUTS:
%   img        -- semnnalul 2D de intrare: imagini alb-negru (NL*NC) SAU imagini color (NL*NC*3)   
%
% OUTPUT:
%   img_THD    -- coeficientii transformatei Hartley NL*NC SAU NL*NC*3 
%   NL, NC     -- dimensiunea semnalului original
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Hartley 2D")


NL=size(img,1);
NC=size(img,2);


numDimensions = ndims(img);

if numDimensions == 2
    % Grayscale image (2D)
    disp("Calcul forma matriceala");
    %formez matricile de transformare pentru forma matriceala
    nl = 0:(NL-1);
    nc = 0:(NC-1);
    kl=nl;
    kc=nc;
    
    [nl,kl] = meshgrid(nl,kl);
    [nc,kc] = meshgrid(nc,kc);
    
    HL= cos(2*pi*nl.*kl/NL) + sin(2*pi*nl.*kl/NL);
    HC= cos(2*pi*nc.*kc/NC) + sin(2*pi*nc.*kc/NC);


    disp("Obtinere coeficienti");
    %obtinerea coeficientilor
    img_double = double(img);
    img_THD = HL*img_double*HC;

elseif numDimensions == 3
    % RGB image (3D)
    
    % Extract each color channel
    img_red = img(:, :, 1);
    img_green = img(:, :, 2);
    img_blue = img(:, :, 3);

    disp("Calcul forma matriceala");
    %formez matricile de transformare pentru forma matriceala
    nl = 0:(NL-1);
    nc = 0:(NC-1);
    kl=nl;
    kc=nc;
    
    [nl,kl] = meshgrid(nl,kl);
    [nc,kc] = meshgrid(nc,kc);
    
    HL= cos(2*pi*nl.*kl/NL) + sin(2*pi*nl.*kl/NL);
    HC= cos(2*pi*nc.*kc/NC) + sin(2*pi*nc.*kc/NC);


    disp("Obtinere coeficienti");
    %obtinerea coeficientilor
    
    img_THD_r = HL*double(img_red)*HC;
    img_THD_g = HL*double(img_green)*HC;
    img_THD_b = HL*double(img_blue)*HC;

    img_THD = cat(3, img_THD_r, img_THD_g, img_THD_b);

else
    error('Unsupported image format. Input must be either a 2D or 3D matrix.');
end


%% SOLUTION END %%

end