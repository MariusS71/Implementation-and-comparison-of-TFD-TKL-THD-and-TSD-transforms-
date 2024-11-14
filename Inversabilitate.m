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
    ylabel('Audio magnitudes');
   

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
   
    chunkSize = 10000;
    signalLength = length(signal);

    % Initialize the result vector
    result = zeros(size(signal));
    signal_tfd = zeros(size(signal));

    % Apply the transform in chunks
    for startIdx = 1:chunkSize:signalLength
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the chunk
        chunk = signal(startIdx:endIdx);
        
        % Transform the chunk
        disp("Aplicare Fourier Discret")
        tfd_chunk = fft(chunk);
        
        % Store the result in the corresponding indices
        result(startIdx:endIdx) = tfd_chunk;
    end

    % Invert the transform in chunks
    for startIdx = 1:chunkSize:signalLength
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the transformed chunk
        tfd_chunk = result(startIdx:endIdx);
        
        % Inverse transform the chunk
        disp("Aplicare Fourier Discret invers")
        signal_tfd_chunk = ifft(tfd_chunk);
        
        % Store the result in the corresponding indices
        signal_tfd(startIdx:endIdx) = signal_tfd_chunk;
    end


    figure('Name','Semnal Fourier', 'Position', [500 100 900 600])
    plot(signal_tfd)
    xlim([0, length(signal_tfd)]);
    ylim([min(signal_tfd)*1.1, max(signal_tfd)*1.1]);
    title('Transformata Fourier Semnal Reconstruit');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');

    figure('Name','Eroare Fourier', 'Position', [500 100 900 600])
    plot(signal-signal_tfd)
    xlim([0, length(signal)]);
    ylim([min(signal-signal_tfd)*1.1, max(signal-signal_tfd)*1.1]);
    title('Transformata Fourier Eroare');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');
    % Calculate the norm value
    errorNorm_TFD = norm(signal - signal_tfd);
    text(signalLength*0.45, max(signal - signal_tfd), ['Norma: ' num2str(errorNorm_TFD)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);

elseif check == 2 && numDimensions == 2

    disp("Aplicare Fourier Discret")
    tfd_result = fft(signal);

    disp("Aplicare Fourier Discret invers")
    tfd_signal = ifft(tfd_result);

    figure('Name','Imagine TFD');
    imshow(uint8(tfd_signal));
    title('Imagine reconstruita Fourier Discret');



    error_tfd=double(signal)-tfd_signal;
    errorNorm_TFD=norm(error_tfd);

    N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare Fourier')
    surf(cols, rows, error_tfd, 'EdgeColor', 'none');
    title('Transformata Fourier Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_tfd(:))*0.45, 0.2, ['Norma: ' num2str(errorNorm_TFD)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);



elseif check == 2 && numDimensions == 3

    disp("Aplicare Fourier Discret")
    tfd_result = fft(signal);

    disp("Aplicare Fourier Discret invers")
    tfd_signal = ifft(tfd_result);

    figure('Name','Imagine TFD');
    imshow(uint8(tfd_signal));
    title('Imagine reconstruita Fourier Discret');


    error_tfd = double(signal) - tfd_signal;
    errorNorm_TFD = norm(error_tfd(:));

    error_tfd = sum(error_tfd, 3);
    
    N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare Fourier')
    surf(cols, rows, error_tfd, 'EdgeColor', 'none');
    title('Transformata Fourier Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_tfd(:))*0.45, 0.2, ['Norma: ' num2str(errorNorm_TFD)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);


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

        % Inverse transform the chunk
        signal_tkl_chunk = inv_TKL_1D(tkl_chunk, Vm, m);
        
        % Store the result in the corresponding indices
        signal_tkl(startIdx:endIdx) = signal_tkl_chunk;
    end


    figure('Name','Semnal TKL', 'Position', [500 100 900 600])
    plot(signal_tkl)
    xlim([0, length(signal_tkl)]);
    ylim([min(signal_tkl)*1.1, max(signal_tkl)*1.1]);
    title('TKL Semnal Reconstruit');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');

    figure('Name','Eroare TKL', 'Position', [500 100 900 600])
    plot(signal-signal_tkl)
    xlim([0, length(signal)]);
    ylim([min(signal-signal_tkl)*1.1, max(signal-signal_tkl)*1.1]);
    title('TKL Eroare');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');
    % Calculate the norm value
    errorNorm_tkl = norm(signal - signal_tkl);
    text(signalLength*0.45, max(signal - signal_tkl), ['Norma: ' num2str(errorNorm_tkl)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);

elseif check == 2 && numDimensions == 2


    [Am, m, Vm] = TKL_2D(signal);

    tkl_signal = inv_TKL_2D(Am, Vm, m);

    figure('Name','Imagine TKL');
    imshow(uint8(tkl_signal));
    title('Imagine reconstruita TKL');



    error_tkl=double(signal)-tkl_signal;
    errorNorm_TKL=norm(error_tkl);
    N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare TKL')
    surf(cols, rows, error_tkl, 'EdgeColor', 'none');
    title('TKL Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_tkl(:))*0.85, -0.2, ['Norma: ' num2str(errorNorm_TKL)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);



elseif check == 2 && numDimensions == 3

    [Am, m, Vm] = TKL_2D(signal);

    tkl_signal = inv_TKL_2D(Am, Vm, m);

    figure('Name','Imagine TKL');
    imshow(uint8(tkl_signal));
    title('Imagine reconstruita TKL');



    error_tkl=double(signal)-tkl_signal;
    errorNorm_TKL=norm(error_tkl(:));
    error_tkl = sum(error_tkl, 3);
    N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare TKL')
    surf(cols, rows, error_tkl, 'EdgeColor', 'none');
    title('TKL Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_tkl(:))*0.85, -2, ['Norma: ' num2str(errorNorm_TKL)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);



end



%% Transformata Hartley

if check == 1
   
    chunkSize = 10000;
    signalLength = length(signal);

    % Initialize the result vector
    signal_thd = zeros(size(signal));

    % Apply the transform in chunks
    i=1;
    for startIdx = 1:chunkSize:signalLength
        disp(['Iteratia ' num2str(i)]);
        i=i+1;
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the chunk
        chunk = signal(startIdx:endIdx);
        
        % Transform the chunk
        [thd_chunk, N] = THD_1D(chunk);

        % Inverse transform the chunk
        signal_thd_chunk = inv_THD_1D(thd_chunk, N);
        
        % Store the result in the corresponding indices
        signal_thd(startIdx:endIdx) = signal_thd_chunk;
    end


    figure('Name','Semnal THD', 'Position', [500 100 900 600])
    plot(signal_thd)
    xlim([0, length(signal_thd)]);
    ylim([min(signal_thd)*1.1, max(signal_thd)*1.1]);
    title('THD Semnal Reconstruit');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');

    figure('Name','Eroare THD', 'Position', [500 100 900 600])
    plot(signal-signal_thd)
    xlim([0, length(signal)]);
    ylim([min(signal-signal_thd)*1.1, max(signal-signal_thd)*1.1]);
    title('THD Eroare');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');
    % Calculate the norm value
    errorNorm_thd = norm(signal - signal_thd);
    text(signalLength*0.45, max(signal - signal_thd), ['Norma: ' num2str(errorNorm_thd)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);

elseif check == 2 && numDimensions == 2


    [img_THD, NL, NC] = THD_2D(signal);

    thd_signal = inv_THD_2D(img_THD, NL, NC);

    figure('Name','Imagine THD');
    imshow(uint8(thd_signal));
    title('Imagine reconstruita THD');



    error_thd=double(signal)-thd_signal;
    errorNorm_THD=norm(error_thd);
       N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare THD')
    surf(cols, rows, error_thd, 'EdgeColor', 'none');
    title('THD Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_thd(:))*0.85, -0.9, ['Norma: ' num2str(errorNorm_THD)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);


elseif check == 2 && numDimensions == 3

    [img_THD, NL, NC] = THD_2D(signal);

    thd_signal = inv_THD_2D(img_THD, NL, NC);

    figure('Name','Imagine THD');
    imshow(uint8(thd_signal));
    title('Imagine reconstruita THD');



    error_thd=double(signal)-thd_signal;
    errorNorm_THD=norm(error_thd(:));
    error_thd = sum(error_thd, 3);

    N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare THD')
    surf(cols, rows, error_thd, 'EdgeColor', 'none');
    title('THD Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_thd(:))*0.85, -0.9, ['Norma: ' num2str(errorNorm_THD)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);


end

%% Transformata Slang

if check == 1
   
    chunkSize = 4096;
    signalLength = length(signal);

    % Initialize the result vector
    nr = signalLength/4096;

    signal_tds = zeros(signalLength, 1);


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
            if(2^nextPowerOf2<4096)
                S = Slant_matrix(nextPowerOf2);  
            end
        end

        i=i+1;
        endIdx = min(startIdx + chunkSize - 1, signalLength);
        
        % Extract the chunk
        chunk = signal(startIdx:endIdx);
        
        % Transform the chunk
        tds_chunk = TDS_1D(chunk, S);

        % Inverse transform the chunk
        signal_tds_chunk = inv_TDS_1D(tds_chunk,S);
        
        if i-1 == ceil(nr)
            endIdx = signalLength;
            signal_tds_chunk = signal_tds_chunk(1:endIdx-startIdx+1);

        end

        % Store the result in the corresponding indices
        signal_tds(startIdx:endIdx) = signal_tds_chunk;
    end

    signal_tds = signal_tds(1:signalLength);

    figure('Name','Semnal TDS', 'Position', [500 100 900 600])
    plot(signal_tds)
    xlim([0, length(signal_tds)]);
    ylim([min(signal_tds)*1.1, max(signal_tds)*1.1]);
    title('TDS Semnal Reconstruit');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');

    figure('Name','Eroare TDS', 'Position', [500 100 900 600])
    plot(signal-signal_tds)
    xlim([0, length(signal)]);
    ylim([min(signal-signal_tds)*1.1, max(signal-signal_tds)*1.1]);
    title('THD Eroare');
    xlabel(['Time (10^{' char('-4') '} seconds)']);
    ylabel('Audio magnitudes');
    % Calculate the norm value
    errorNorm_tds = norm(signal - signal_tds);
    text(signalLength*0.45, max(signal - signal_tds), ['Norma: ' num2str(errorNorm_tds)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);

elseif check == 2 && numDimensions == 2


    nextPowerOf2 = nextpow2(max(size(signal,1), size(signal,2)));
    S = Slant_matrix(nextPowerOf2);
    
    coef = TDS_2D(double(signal), S);
    
    
    [rows,cols] = size(signal,1,2);
    
    inv = inv_TDS_2D(coef,S, rows, cols);

    figure('Name','Imagine TDS');
    imshow(uint8(inv));
    title('Imagine reconstruita TDS');



    error_tds=double(signal)-inv;
    errorNorm_TDS=norm(error_tds);
        N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare TDS')
    surf(cols, rows, error_tds, 'EdgeColor', 'none');
    title('TDS Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_tds(:))*0.85, -0.2, ['Norma: ' num2str(errorNorm_TDS)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);



elseif check == 2 && numDimensions == 3

    nextPowerOf2 = nextpow2(max(size(signal,1), size(signal,2)));
    S = Slant_matrix(nextPowerOf2);
    
    coef = TDS_2D(double(signal), S);
    
    
    [rows,cols] = size(signal,1,2);
    
    inv = inv_TDS_2D(coef,S, rows, cols);

    figure('Name','Imagine TDS');
    imshow(uint8(inv));
    title('Imagine reconstruita TDS');



    error_tds=double(signal)-inv;
    errorNorm_TDS=norm(error_tds(:));
    error_tds = sum(error_tds, 3);

    N = size(signal);

    
    % Create the frequency axes for rows and columns
    frequencyRows = linspace(0, N(1), N(1)) ;
    frequencyCols = linspace(0, N(2), N(2)) ;
    
    [rows, cols] = meshgrid(frequencyCols, frequencyRows);

    figure('Name','Eroare TDS')
    surf(cols, rows, error_tds, 'EdgeColor', 'none');
    title('TDS Eroare');
    xlabel('Columns');
    ylabel('Rows');
    zlabel("Eroare pixel")
    text(max(error_tds(:))*0.85, -0.2, ['Norma: ' num2str(errorNorm_TDS)], ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);



end
