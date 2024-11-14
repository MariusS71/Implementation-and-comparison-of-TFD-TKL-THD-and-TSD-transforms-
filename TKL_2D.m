function [Am, m, Vm] = TKL_2D(img)
%% Transformata Karhunen-Loeve 2D
% INPUTS:
%   img         -- imaginea originala
%
% OUTPUT:
%   Am          -- coeficientii transformarii TKL
%   Vm          -- vectorii proprii
%   m           -- media semnalului 
%
% OBSERVATII:
%  

%% SOLUTION START %%

disp("Transformata Karhunen-Loeve 2D")


numDimensions = ndims(img);

if numDimensions == 2
    % Grayscale image (2D)
    %centrare
    m=mean(img, 'all');
    A = double(img)-m;
    %A = double(img) - mean(img);
    
    disp("Calcul matricei de auto-corelatie")
    % matricea de auto corelatie
    R=A'*A;
    
    % Calcularea valorilor și vectorilor proprii
    disp("Calcul valori și vectori proprii")
    [V, D] = eig(R);
    
    % Sortarea in ordine descrescatoare
    [~, order] = sort(diag(D), 'descend');
    %D = D(order, order);
    V = V(:, order);

   
    %Vm=V;
    Vm=V(:,1:end);
    
    disp("Calcul coeficienti")
    Am=Vm'*A';


elseif numDimensions == 3
    % RGB image (3D)
    
    % Extract each color channel
    img_red = img(:, :, 1);
    img_green = img(:, :, 2);
    img_blue = img(:, :, 3);


    %centrare
    m=[mean(img_red, 'all'), mean(img_green, 'all'), mean(img_blue, 'all')];
    A_red = double(img_red)-m(1);
    A_green = double(img_green)-m(2);
    A_blue = double(img_blue)-m(3);
    
    %A = double(img) - mean(img);
    
    % matricea de auto corelatie
    disp("Calcul matricei de auto-corelatie")
    R_r=A_red'*A_red;
    R_g=A_green'*A_green;
    R_b=A_blue'*A_blue;
    
    % Calcularea valorilor și vectorilor proprii
    disp("Calcul valori și vectori proprii")
    [V_r, D_r] = eig(R_r);
    [V_g, D_g] = eig(R_g);
    [V_b, D_b] = eig(R_b);

    [~, order_r] = sort(diag(D_r), 'descend');
    [~, order_g] = sort(diag(D_g), 'descend');
    [~, order_b] = sort(diag(D_b), 'descend');
    
    
    V_r = V_r(:, order_r);
    V_g = V_g(:, order_g);
    V_b = V_b(:, order_b);
    
    
    Vm_r=V_r(:,1:end);
    Vm_g=V_g(:,1:end);
    Vm_b=V_b(:,1:end);
    Vm = cat(3, Vm_r, Vm_g, Vm_b);

    disp("Calcul coeficienti")
    Am_r=Vm_r'*A_red';
    Am_g=Vm_g'*A_green';
    Am_b=Vm_b'*A_blue';
    Am = cat(3, Am_r, Am_g, Am_b);


else
    error('Unsupported image format. Input must be either a 2D or 3D matrix.');
end


%% SOLUTION END %%

end