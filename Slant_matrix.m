function [M] = Slant_matrix(num)
%% Construirea matricii Slant
% INPUTS:
%   num    --  ordinul matricii Slant
%   
%
% OUTPUT:
%   M      -- matricea Slant (dimensiune 2^num x 2^num) ordonata
%             crescator in functie de benzile dominante de frecventa
%
% OBSERVATII:
%  

delta = sqrt(2);

l=2;
delta_new = (delta *  sqrt(2)) / sqrt( 2^(3*(l-1)) * delta^2 + 4);


S = (1/sqrt(2)) * [1 1 ; 1 -1];

S_new = (1/2)* [ 1 1 1 1; 3/sqrt(5)  1/sqrt(5) -1/sqrt(5) -3/sqrt(5);
                1 -1 -1 1; 1/sqrt(5) -3/sqrt(5) 3/sqrt(5) -1/sqrt(5)
                  ];
a=2/sqrt(5);
b=1/sqrt(5);

A = [ 1 0; a b];
B = [ 1 0; -a b];
C = [0 1; -b a];
D = [0 -1; b a];



if (num == 1)  
    M = S;
elseif (num == 2)
    M=S_new;
else 
    l=l+1;

    disp("Calcul S");

    while(l<=num)
       delta = delta_new;
       aux = 2^(3*(l-1)) * delta^2 + 4;
       a=(2^(3*(l-1)/2) * delta) / sqrt(aux);

       b=2/sqrt(aux);

       delta_new = (delta *  sqrt(2)) / sqrt(aux);


        
        S=S_new;

        A = [ A, zeros(size(A,1)); zeros(size(A,1)), eye(size(A,1))];
        B = [ B, zeros(size(B,1)); zeros(size(B,1)), eye(size(B,1))];
        C = [ C, zeros(size(C,1)); zeros(size(C,1)), eye(size(C,1))];
        D = [ D, zeros(size(D,1)); zeros(size(D,1)), -eye(size(D,1))];

        A(2,1)=a;
        A(2,2)=b;

        B(2,1)=-a;
        B(2,2)=b;

        C(2,1)=-b;
        C(2,2)=a;

        D(2,1)=b;
        D(2,2)=a;


        S_new= (1/sqrt(2)) * [ A, B; C, D] * [S, zeros(size(S,1));
                                              zeros(size(S,1)), S ];
        

        l=l+1;
    end

    S=S_new; % Matricea Slant

    % ordonarea crescatoare a liniilor in functie de banda dominanta de
    % frecventa
    
    disp("Ordonare S");

    % Parameters
    Fs = 1000; % Sampling frequency
    
    % Initialize arrays to store dominant frequencies and corresponding times
    dominant_frequencies = zeros(size(S, 1), 1);
    
    % Calculate spectrum for each line
    for i = 1:size(S, 1)
        signal = S(i, :);
    
        N = 100 * length(signal); % Increase the number of points
        % Calculate FFT
    
    
        signal_fft= fft(signal , N);
        amplitude_spectrum = abs(signal_fft/N);
        amplitude_spectrum = amplitude_spectrum(1:N/2+1);
        frequencies = Fs * (0:(N/2))/N;
    
        % Plot the amplitude spectrum
        %figure
        %plot(frequencies, 2*amplitude_spectrum);
    
        % Find the dominant frequency
        [~, max_idx] = max(amplitude_spectrum);
        dominant_frequencies(i) = frequencies(max_idx);
    
    end
    
    % Order signals based on dominant time
    [~, order] = sort(dominant_frequencies);
    
    % Reorder the matrix of signals
    M = S(order, :);


end


end