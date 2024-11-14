%% Program de testare a inversabilitatii transformatelor (TFD, TKL, THD, TSD)
% !! Se recomanda rularea programului pe sectiuni. Prima sectiune este
% destinata citirii semnalului. Se recomanda rularea primei sectiuni
% inainte de fiecare transformata pentru ca numele variabilelor sa nu se
% suprapuna !!


clear; clc; close all;

check = 0;  % Initialize check variable

while check == 0
    % Prompt the user for the file name or path
    fileName = input('Enter the name or path of the signal file: ', 's');

    % Check if the file exists
    if exist(fileName, 'file') == 2
        % Determine the file type based on the file extension
        [~, ~, fileExtension] = fileparts(fileName);

        % Load the signal based on the file type
        switch lower(fileExtension)
            case {'.wav', '.mp3', '.flac', '.ogg'} % Audio file types
                signal = audioread(fileName);
                disp(['Audio file "' fileName '" loaded successfully.']);
                check = 1;

            case {'.jpg', '.png', '.bmp', '.gif'} % Image file types
                signal = imread(fileName);
                disp(['Image file "' fileName '" loaded successfully.']);
                check = 2;
                numDimensions = ndims(signal);

            otherwise
                disp('Unsupported file type. Please provide a valid audio or image file.');
        end
    else
        disp(['File "' fileName '" not found. Please enter a valid file name.']);
    end
end


if check == 1
    figure('Name','Semnal Original', 'Position', [500 100 900 600])
    plot(signal)
    xlim([0, length(signal)]);
    ylim([min(signal)*1.1, max(signal)*1.1]);
    title('Semnal Original');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Amplitude');
   

elseif check == 2 && numDimensions == 2
    figure('Name','Semnal Original');
    imshow(signal);
    title('Imagine originala');

elseif check == 2 && numDimensions == 3
    figure('Name','Semnal Original');
    imshow(signal);
    title('Imagine originala');

end    

%% Transformata Fourier Discreta

if check == 1
   
    
    Fs = 1;  
    
    % Apply FFT to the signal
    tfd = fft(signal);
    
    % Calculate the frequency axis
    N = length(tfd);
    frequencies = linspace(0, 0.5, N/2 + 1) * Fs;
    
    % Calculate the magnitude spectrum in decibels
    magnitude_dB = 20 * log10(abs(tfd(1:N/2 + 1)));
    
    % Plot the Spectral Power in dB using Normalized Frequency
    figure;
    plot(frequencies, magnitude_dB);
    title('Spectrum of signal');
    xlabel('Normalized Frequency');
    ylabel('Spectral Power (dB)');
    
    
    
    % Calculate the total energy of the signal
    totalEnergy = sum(10.^(magnitude_dB/20));
    
    % Sort magnitude_dB in descending order
    sortedMagnitude_dB = sort(magnitude_dB, 'descend');
    %sortedMagnitude_dB = magnitude_dB;
    
    % Calculate the desired percentage of total energy (e.g., 45%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < N && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + 10^(sortedMagnitude_dB(numElements)/20);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude_dB)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);

elseif check == 2 && numDimensions == 2

    Fs = 1;  
    
    tfd = fft2(signal);

    % Calculate the magnitude spectrum in decibels
    N = size(tfd);
    magnitude_dB = 20 * log10(abs(tfd(1:N(1)/2 + 1, 1:N(2)/2 + 1)));
    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 0.5, N(1)/2 + 1) * Fs;
    frequencyCols = linspace(0, 0.5, N(2)/2 + 1) * Fs;
    
    % Create a meshgrid for visualization (optional)
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power in dB using Normalized Frequency
    figure;
    surf(cols, rows, magnitude_dB, 'EdgeColor', 'none');
    title('Spectrum of Image');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power (dB)');
    
    
    
   % Calculate the total energy of the image
    totalEnergy = sum(10.^(magnitude_dB(:)/20));
    
    % Sort magnitude_dB in descending order
    sortedMagnitude_dB = sort(magnitude_dB(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude_dB) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + 10^(sortedMagnitude_dB(numElements)/20);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude_dB)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
    percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);





elseif check == 2 && numDimensions == 3
    
    Fs=1;

    tfd_R = fft2(signal(:, :, 1));  % Red channel
    tfd_G = fft2(signal(:, :, 2));  % Green channel
    tfd_B = fft2(signal(:, :, 3));  % Blue channel
    
    N = size(signal);
    
    % Calculate the magnitude spectrum in decibels for each channel
    magnitude_dB_R = 20 * log10(abs(tfd_R(1:N(1)/2 + 1, 1:N(2)/2 + 1)));
    magnitude_dB_G = 20 * log10(abs(tfd_G(1:N(1)/2 + 1, 1:N(2)/2 + 1)));
    magnitude_dB_B = 20 * log10(abs(tfd_B(1:N(1)/2 + 1, 1:N(2)/2 + 1)));

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 0.5, N(1)/2 + 1) * Fs;
    frequencyCols = linspace(0, 0.5, N(2)/2 + 1) * Fs;
    
    % Create meshgrids for visualization (optional)
    [rows_R, cols_R] = meshgrid(frequencyCols, frequencyRows);
    [rows_G, cols_G] = meshgrid(frequencyCols, frequencyRows);
    [rows_B, cols_B] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power in dB for each channel using Normalized Frequency
    
    figure
    surf(cols_R, rows_R, magnitude_dB_R, 'EdgeColor', 'none');
    title('Red Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power (dB)');
    
     figure
    surf(cols_G, rows_G, magnitude_dB_G, 'EdgeColor', 'none');
    title('Green Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power (dB)');
    
     figure
    surf(cols_B, rows_B, magnitude_dB_B, 'EdgeColor', 'none');
    title('Blue Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power (dB)');

    
    % Combine the magnitudes from each channel (e.g., using average)
    magnitude_dB = (magnitude_dB_R + magnitude_dB_G + magnitude_dB_B)/3 ;
    
    % Calculate the total energy of the image
    totalEnergy = sum(10.^(magnitude_dB(:)/20));
    
    % Sort magnitude_dB in descending order
    sortedMagnitude_dB = sort(magnitude_dB(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude_dB) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + 10^(sortedMagnitude_dB(numElements)/20);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude_dB)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);

    

end


%% Transformata Karhunen-Loeve



if check == 1
    
    
    chunkSize = 10000;
    signalLength = length(signal);

    % Initialize the result vector
    signal_tkl = zeros(size(signal));

    % Apply the transform in chunks
    i=1;
    for startIdx = 1:chunkSize:signalLength
        disp(['Iteratia ' num2str(i)]);
        i=i+1;
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the chunk
        chunk = signal(startIdx:endIdx);
        
        % Transform the chunk
        [tkl_chunk, m, Vm] = TKL_1D(chunk);      
        
        % Store the result in the corresponding indices
        signal_tkl(startIdx:endIdx) = tkl_chunk;
    end


    % Calculate the frequency axis
    Fs=1;
    N = length(signal_tkl);
    
    magnitude = abs(signal_tkl);
    
    % Calculate the total energy of the signal
    totalEnergy = sum(magnitude);


    frequencies = linspace(0, 1, N) * Fs;
    
    
    % Plot the Spectral Power in dB using Normalized Frequency
    figure;
    plot(frequencies, signal_tkl);
    title('Spectrum of signal');
    xlabel('Normalized Frequency');
    ylabel('Spectral Power');
    
    
    
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude, 'descend');
    

    desiredPercentage = 99;

    % Calculate the desired percentage of total energy (e.g., 45%)
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        %cumulativeEnergy = cumulativeEnergy + 10^(sortedMagnitude_dB(numElements)/20);
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / N) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);





elseif check == 2 && numDimensions == 2
     
    

    [Am, m, Vm] = TKL_2D(signal);

    % Calculate the frequency axis
    Fs = 1;  
    
    N = size(Am);

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 1, N(1)) * Fs;
    frequencyCols = linspace(0, 1, N(2)) * Fs;
    
    % Create a meshgrid for visualization (optional)
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power using Normalized Frequency
    figure;
    surf(cols, rows, abs(Am), 'EdgeColor', 'none');
    title('Spectrum of Image');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
    
    
   % Calculate the total energy of the image
    magnitude = abs(Am);

    totalEnergy = sum(magnitude(:));
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while (numElements < length(sortedMagnitude)) && (cumulativeEnergy < desiredEnergy)
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);





elseif check == 2 && numDimensions == 3


    [Am, m, Vm] = TKL_2D(signal);

    % Calculate the frequency axis
    Fs=1;

    Am_R = Am(:, :, 1);  % Red channel
    Am_G = Am(:, :, 2);  % Green channel
    Am_B = Am(:, :, 3);  % Blue channel
    
    N = size(signal);
    

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 1, N(2)) * Fs;
    frequencyCols = linspace(0, 1, N(1)) * Fs;
    
    % Create meshgrids for visualization (optional)
    [rows_R, cols_R] = meshgrid(frequencyCols, frequencyRows);
    [rows_G, cols_G] = meshgrid(frequencyCols, frequencyRows);
    [rows_B, cols_B] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power in dB for each channel using Normalized Frequency
    
    figure
    surf(cols_R, rows_R, Am_R, 'EdgeColor', 'none');
    title('Red Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
     figure
    surf(cols_G, rows_G, Am_G, 'EdgeColor', 'none');
    title('Green Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
     figure
    surf(cols_B, rows_B, Am_B, 'EdgeColor', 'none');
    title('Blue Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');

    % Calculate the magnitude spectrum for each channel
    magnitude_R = abs(Am_R);
    magnitude_G = abs(Am_G);
    magnitude_B = abs(Am_B);
    
    % Combine the magnitudes from each channel (e.g., using average)
    magnitude = (magnitude_R + magnitude_G + magnitude_B)/3 ;
    
    % Calculate the total energy of the image
    totalEnergy = sum(magnitude(:));
    
    % Sort magnitude_dB in descending order
    sortedMagnitude = sort(magnitude(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);
    disp(numElements)


end



%% Transformata Hartley

if check == 1
   
    chunkSize = 20000;
    signalLength = length(signal);

    % Initialize the result vector
    thd = zeros(size(signal));

    % Apply the transform in chunks
    i=1;
    for startIdx = 1:chunkSize:signalLength
        disp(['Iteratia ' num2str(i)]);
        i=i+1;
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the chunk
        chunk = signal(startIdx:endIdx);
        
        % Transform the chunk
        [thd_chunk, ~] = THD_1D(chunk);
        
        % Store the result in the corresponding indices
        thd(startIdx:endIdx) = thd_chunk;
    end


   Fs = 1;  % Replace with your actual sampling rate
    

    % Calculate the frequency axis
    N = length(thd);
    frequencies = linspace(0, 1, N);
    

    
    % Plot the Spectral Power using Normalized Frequency
    figure;
    plot(frequencies, thd);
    title('Spectrum of signal (Hartley Transform)');
    xlabel('Normalized Frequency');
    ylabel('Spectral Power');

    magnitude = abs(thd);
    
    % Calculate the total energy of the signal
    totalEnergy = sum(magnitude);
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude, 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 45;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / N) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
    percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);

elseif check == 2 && numDimensions == 2


    [img_THD, NL, NC] = THD_2D(signal);

    Fs = 1;  
    
    N = size(img_THD);

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 1, N(1)) * Fs;
    frequencyCols = linspace(0, 1, N(2)) * Fs;
    
    % Create a meshgrid for visualization (optional)
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power using Normalized Frequency
    figure;
    surf(cols, rows, img_THD, 'EdgeColor', 'none');
    title('Spectrum of Image');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
    
    
   % Calculate the total energy of the image
    magnitude = abs(img_THD);

    totalEnergy = sum(magnitude(:));
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);



elseif check == 2 && numDimensions == 3

    [img_THD, NL, NC] = THD_2D(signal);

    Fs=1;

    img_THD_R = img_THD(:, :, 1);  % Red channel
    img_THD_G = img_THD(:, :, 2);  % Green channel
    img_THD_B = img_THD(:, :, 3);  % Blue channel
    
    N = size(signal);
    

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 1, N(1)) * Fs;
    frequencyCols = linspace(0, 1, N(2)) * Fs;
    
    % Create meshgrids for visualization (optional)
    [rows_R, cols_R] = meshgrid(frequencyCols, frequencyRows);
    [rows_G, cols_G] = meshgrid(frequencyCols, frequencyRows);
    [rows_B, cols_B] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power in dB for each channel using Normalized Frequency
    
    figure
    surf(cols_R, rows_R, img_THD_R, 'EdgeColor', 'none');
    title('Red Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
     figure
    surf(cols_G, rows_G, img_THD_G, 'EdgeColor', 'none');
    title('Green Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
     figure
    surf(cols_B, rows_B, img_THD_B, 'EdgeColor', 'none');
    title('Blue Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');

    % Calculate the magnitude spectrum for each channel
    magnitude_R = abs(img_THD_R);
    magnitude_G = abs(img_THD_G);
    magnitude_B = abs(img_THD_B);
    
    % Combine the magnitudes from each channel (e.g., using average)
    magnitude = (magnitude_R + magnitude_G + magnitude_B)/3 ;
    
    % Calculate the total energy of the image
    totalEnergy = sum(magnitude(:));
    
    % Sort magnitude_dB in descending order
    sortedMagnitude = sort(magnitude(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);

end

%% Transformata Slang

if check == 1
   
    chunkSize = 4096;
    signalLength = length(signal);

    % Initialize the result vector
    nr = signalLength/4096;
    if nr > 1
        tds = zeros(4096 * floor(nr) + 2^nextpow2(signalLength - 4096 * floor(nr)), 1);
    else
        tds = zeros(2^nextpow2(signalLength), 1);
    end

    if signalLength >= 4096
        nextPowerOf2 = nextpow2(4096);
        S = Slant_matrix(nextPowerOf2);

    else
        nextPowerOf2 = nextpow2(signalLength);
        S = Slant_matrix(nextPowerOf2);
    end


    % Apply the transform in chunks
     
    i=1;
    for startIdx = 1:chunkSize:signalLength
        disp(['Iteratia ' num2str(i)]);

        if i == ceil(nr) %% daca suntem in ultima iteratiei a buclei dimensiune chunk <= 4096
            nextPowerOf2 = nextpow2(signalLength - endIdx);
            if 2^nextPowerOf2<4096
                S = Slant_matrix(nextPowerOf2);  
            end
        end

        i=i+1;
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the chunk
        chunk = signal(startIdx:endIdx);
        
        % Transform the chunk
        tds_chunk = TDS_1D(chunk, S);

        
        if i-1 == ceil(nr)
            endIdx = 4096 * floor(nr) + 2^nextpow2(signalLength - 4096 * floor(nr));
        end

        % Store the result in the corresponding indices
        tds(startIdx:endIdx) = tds_chunk;
    end

   Fs = 1;  % Replace with your actual sampling rate
    

    % Calculate the frequency axis
    N = length(tds);
    frequencies = linspace(0, 1, N);
    

    
    % Plot the Spectral Power using Normalized Frequency
    figure;
    plot(frequencies, tds);
    title('Spectrum of signal (Slant Transform)');
    xlabel('Normalized Frequency');
    ylabel('Spectral Power');

    magnitude = abs(tds);
    
    % Calculate the total energy of the signal
    totalEnergy = sum(magnitude);
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude, 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < N && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / N) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
    percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);

elseif check == 2 && numDimensions == 2


    nextPowerOf2 = nextpow2(max(size(signal,1), size(signal,2)));
    S = Slant_matrix(nextPowerOf2);
    
    img_TDS = TDS_2D(double(signal), S);


    Fs = 1;  
    
    N = size(img_TDS);

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 1, N(1)) * Fs;
    frequencyCols = linspace(0, 1, N(2)) * Fs;
    
    % Create a meshgrid for visualization (optional)
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power using Normalized Frequency
    figure;
    surf(cols, rows, img_TDS, 'EdgeColor', 'none');
    title('Spectrum of Image');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
    
    
   % Calculate the total energy of the image
    magnitude = abs(img_TDS);

    totalEnergy = sum(magnitude(:));
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);



elseif check == 2 && numDimensions == 3

    nextPowerOf2 = nextpow2(max(size(signal,1), size(signal,2)));
    S = Slant_matrix(nextPowerOf2);
    
    img_TDS = TDS_2D(double(signal), S);

    Fs=1;

    img_TDS_R = img_TDS(:, :, 1);  % Red channel
    img_TDS_G = img_TDS(:, :, 2);  % Green channel
    img_TDS_B = img_TDS(:, :, 3);  % Blue channel
    
    N = size(img_TDS_R);
    

    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, 1, N(1)) * Fs;
    frequencyCols = linspace(0, 1, N(2)) * Fs;
    
    % Create meshgrids for visualization (optional)
    [rows_R, cols_R] = meshgrid(frequencyCols, frequencyRows);
    [rows_G, cols_G] = meshgrid(frequencyCols, frequencyRows);
    [rows_B, cols_B] = meshgrid(frequencyCols, frequencyRows);
    
    % Plot the Spectral Power in dB for each channel using Normalized Frequency
    
    figure
    surf(cols_R, rows_R, img_TDS_R, 'EdgeColor', 'none');
    title('Red Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
     figure
    surf(cols_G, rows_G, img_TDS_G, 'EdgeColor', 'none');
    title('Green Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');
    
     figure
    surf(cols_B, rows_B, img_TDS_B, 'EdgeColor', 'none');
    title('Blue Channel Spectrum');
    xlabel('Normalized Frequency (Columns)');
    ylabel('Normalized Frequency (Rows)');
    zlabel('Spectral Power');

    % Calculate the magnitude spectrum for each channel
    magnitude_R = abs(img_TDS_R);
    magnitude_G = abs(img_TDS_G);
    magnitude_B = abs(img_TDS_B);
    
    % Combine the magnitudes from each channel (e.g., using average)
    magnitude = (magnitude_R + magnitude_G + magnitude_B)/3 ;
    
    % Calculate the total energy of the image
    totalEnergy = sum(magnitude(:));
    
    % Sort magnitude in descending order
    sortedMagnitude = sort(magnitude(:), 'descend');
    
    % Calculate the desired percentage of total energy (e.g., 55%)
    desiredPercentage = 99;
    desiredEnergy = totalEnergy * desiredPercentage / 100;
    
    % Calculate the number of elements needed to reach the desired energy
    cumulativeEnergy = 0;
    numElements = 0;
    while numElements < length(sortedMagnitude) && cumulativeEnergy < desiredEnergy
        numElements = numElements + 1;
        cumulativeEnergy = cumulativeEnergy + sortedMagnitude(numElements);
    end
    
    % Calculate the percentage of the first elements that encapsulate the desired energy
    percentageFirstElements = (numElements / numel(magnitude)) * 100;
    
    fprintf('Percentage of the first %.2f%% elements that encapsulate %.2f%% of the energy: %.2f%%\n', ...
        percentageFirstElements, desiredPercentage, cumulativeEnergy / totalEnergy * 100);

end

