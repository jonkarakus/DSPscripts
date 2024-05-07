%This MATLAB script, authored by Jon Karakus, implements a speech coding system 
% based on Linear Predictive Coding (LPC) and quantization techniques. The code 
% aims to compress speech signals by exploiting their predictability and reducing 
% the number of bits required for representation.

%-----------Speech Signal Quantization code-------------------------

% Initialize variables
frameSize  = 800;
fftLen     = 2048;
analysis   = true;

% Array of quantization bits to test
quantization_bits = [2, 4, 6, 8];

audioReader = dsp.AudioFileReader('speech_file.wav', 'SamplesPerFrame', frameSize, ...
            'OutputDataType', 'double');

% Process and plot for each quantization bit depth
for bits = quantization_bits
    % Reset audio reader to start from the beginning of the file
    reset(audioReader);
    
    % Create figure for plots
    figure('Name', ['Quantization with ' num2str(bits) ' bits'], 'NumberTitle', 'off');
    tiledlayout(2, 1);  % Two subplots: original and quantized signal
    
    % Stream Processing Loop
    while ~isDone(audioReader)
        % Read audio input
        sig = audioReader();

        % Analysis
        sigpreem = preEmphasisFilter(sig);  % Pre-emphasize the signal
        write(signalBuffer, sigpreem);  % Write the pre-emphasized signal into the buffer
        sigbuf = read(signalBuffer, 2 * frameSize, frameSize);  % Read from the buffer with overlap
        hammingwin = hamming(2 * frameSize);  % Create a Hamming window
        sigwin = hammingwin .* sigbuf;  % Apply the window to the buffered signal
        sigacf = xcorr(sigwin, 12, 'biased');  % Compute the autocorrelation
        sigacf = sigacf(13:end);  % Select the right part of the autocorrelation
        [sigA, ~, sigK] = levinson(sigacf);  % Levinson-Durbin recursion to compute LPC
        siglpc = analysisFilter(sigpreem, sigK);  % Filter the signal with the LPC to get the prediction error

        % Scale the prediction error to avoid saturation and underflow
        maxVal = max(abs(siglpc));
        if maxVal == 0
            maxVal = 1;  % Avoid division by zero if siglpc is a zero signal
        end
        siglpc_scaled = siglpc / maxVal;

        % Quantize the prediction error
        sigq = fixq(siglpc_scaled, bits, 'round', 'sat');

        % Rescale the quantized signal back to the original range
        sigq_rescaled = sigq * maxVal;

        % Synthesis
        synthesisFilter.ReflectionCoefficients = sigK.';
        sigsyn = synthesisFilter(sigq);          
        sigout = deEmphasisFilter(sigsyn);   
    
        audioWriter(sigout); % 
        % Play output audio
        audioWriter(sigout);
    end
    
    % Plot the original and quantized signals
    nexttile;
    plot(siglpc);
    title('Original Prediction Error');
    xlim([0 300])

    nexttile;
    plot(sigq_rescaled);
    title(['Quantized Prediction Error with ' num2str(bits) ' bits']);
    xlim([0 300])
    
    % Update scope plots if needed
    % ... [Update the dsp.SpectrumAnalyzer plots if necessary] ...
    
    % Hold on to the figure before moving to the next quantization bit level
    pause(1);  % Pause to update the plots
end

% Release System objects
release(audioReader);
release(audioWriter);
release(scope);
