% Code authored by Jon Karakus to investigate the Least Mean Squares (LMS) adaptation algorithm
% with White Noise and Speech Inputs. This script analyzes the algorithm's performance and
% efficacy in different acoustic scenarios.

%----------------------------Plant Simulation --------------------------

% load('elliptic.mat'); % Loads the SOS matrix and G vector
% [h, t] = impz(SOS, 160, 8000); % Computes the impulse response for 160 samples (20 ms at 8000 Hz)
% h = prod(G) * h; % Scales the impulse response by the product of the G vector
% plot(t, h); % Plots the scaled and truncated impulse response
% title('Scaled Impulse Response Truncated to 20 ms');
% xlabel('Time (samples)');
% ylabel('Amplitude');
% impz(SOS, 160, 8000); % Computes and plots the impulse response for 160 samples without scaling
% title('Impulse Response without Scaling');
% freqz(SOS); % Plots the frequency response of the filter
% title('Frequency Response with Passband Raised by 6.36 dB');

%----------------------LMS Adaptation --------------------------------

M = 21; % Number of coefficients in the adaptive filter
N = 81; 
mu = 0.01; % Step size for the LMS algorithm
[x, Fs] = audioread('speech.wav'); %replce with desired input
numIterations = length(x) - N + 1; % Number of iterations

w = zeros(N, 1);

e = zeros(numIterations, 1);
y = zeros(numIterations, 1);


for k = 1:numIterations
    % Ensure the index does not exceed the bounds of x
    xk = x(k:min(k+N-1, length(x))); % Extract input sample window
    
    % Compute the current output of the adaptive filter
    yk = w' * flipud(xk);
    
    ek = x(k) - yk; % Error computation
    
    % Update the filter coefficients
    w = w + mu * ek * flipud(xk);
    
    % Store the error and output for performance monitoring
    e(k) = ek;
    y(k) = yk;
end


% Plot the error signal
figure;
plot(e);
title('Error Signal');
xlabel('Iteration (k)');
ylabel('Error Magnitude');

% Plot the adaptive filter coefficients
figure;
plot(w);
title('Adaptive Filter Coefficients');
xlabel('Coefficient Index');
ylabel('Coefficient Value');



