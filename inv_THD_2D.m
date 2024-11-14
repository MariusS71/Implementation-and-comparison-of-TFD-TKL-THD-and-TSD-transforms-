function y_img = inv_THD_2D(img_THD, NL, NC)
%% Transformata Hartley 2D inversa
% INPUTS:
%   img_THD    -- coeficientii transformatei Hartley NL*NC SAU NL*NC*3 
%   NL, NC     -- dimensiunea semnalului original
%
% OUTPUT:
%   img        -- semnnalul 2D reconstruit: imagini alb-negru (NL*NC) SAU imagini color (NL*NC*3)   
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Hartley inversa 2D")

numDimensions = ndims(img_THD);

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

    %reconstruire imagine originala
    disp("Obtinere inversa")
    y_img =  (HL * img_THD *HC) / (NL*NC);

   

elseif numDimensions == 3
    % RGB image (3D)
    
    % Extract each color channel
    img_THD_r = img_THD(:, :, 1);
    img_THD_g = img_THD(:, :, 2);
    img_THD_b = img_THD(:, :, 3);

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

    %reconstruire imagine originala
    disp("Obtinere inversa")
    y_img_red =  (HL * img_THD_r *HC) / (NL*NC);
    y_img_green =  (HL * img_THD_g*HC) / (NL*NC);
    y_img_blue =  (HL * img_THD_b*HC) / (NL*NC);
    
    y_img = cat(3, y_img_red, y_img_green, y_img_blue);
    

else
    error('Unsupported image format. Input must be either a 2D or 3D matrix.');
end


%% SOLUTION END %%

end