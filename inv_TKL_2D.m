function B = inv_TKL_2D(Am, Vm, m)
%% Transformata Karhunen-Loeve 2D inversa
% INPUTS:
%  Am           -- coeficientii transformarii TKL
%  Vm           -- vectorii proprii
%  m            -- media semnalului
%
% OUTPUT:
%   B           -- imaginea reconstruita
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Karhunen-Loeve 2D inversa")

numDimensions = ndims(Am);

if numDimensions == 2
    % Grayscale image (2D)
    disp("Reconstruire imagine")
    B = Vm*Am + m; %=~ s
    B=B';
    
    %B=uint8(B);

   

elseif numDimensions == 3

    disp("Reconstruire imagine")
    % Extract each color channel
    Vm_r = Vm(:, :, 1);
    Vm_g = Vm(:, :, 2);
    Vm_b = Vm(:, :, 3);

    Am_r = Am(:, :, 1);
    Am_g = Am(:, :, 2);
    Am_b = Am(:, :, 3);


    % RGB image (3D)
    B_r = Vm_r*Am_r + m(1); 
    B_g = Vm_g*Am_g + m(2); 
    B_b = Vm_b*Am_b + m(3); 
    
    B_r=B_r';
    B_g=B_g';
    B_b=B_b';
    
    B=cat(3, B_r, B_g, B_b);
    %B=uint8(B);
    
    

else
    error('Unsupported image format. Input must be either a 2D or 3D matrix.');
end


%% SOLUTION END %%

end