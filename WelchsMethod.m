
% MATLAB Script for Signal Processing Using Welch's Method
% Author: Jon Karakus
%
% This script is designed to demonstrate the application of the Least Mean Squares (LMS) 
% algorithm to a signal containing white noise and a sinusoidal interference. The script 
% generates a primary signal by superimposing white noise on a sinusoidal noise source.
% It then applies the LMS adaptive filter to estimate the system parameters and 
% subsequently uses Welch's method to estimate the power spectral density (PSD) of the 
% output error signal. This process is repeated for different values of the adaptation 
% step size, mu, to analyze the effect on filter performance.

%-------------------Welch's Method to estimate PSD-----------------
fs = 16000; 
duration = 5; 
amp = 0.3; 
f_noise = 1000; 
N = 16; 
mu = 0.01; 
t = 0:1/fs:duration-1/fs; % Time vector

% Generate the white noise content-carrying signal with the same power as the sine wave
content_signal = normrnd(0, amp/sqrt(2), size(t)); 

noise_source = amp * sin(2 * pi * f_noise * t);

primary_signal = content_signal + noise_source;

reference_signal = noise_source;

% Apply the LMS filter
[a, e, ~] = my_lms(reference_signal, primary_signal, N, mu);

% Estimate the PSD of the output error signal using Welch's method
[pxx, f] = pwelch(e, hamming(256), 128, 256, fs);


mu_values = [0.01, 0.05, 0.1, 0.15, 0.2]; % Example values for mu
for i = 1:length(mu_values)
    mu = mu_values(i);
    [a, e, ~] = my_lms(reference_signal, primary_signal, N, mu);
    [pxx, ~] = pwelch(e, hamming(256), 128, 256, fs);
    figure;
    plot(f, 10*log10(pxx));
    title(['Welchs method PSD of the Error Signal for mu = ', num2str(mu)]);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
end

%-------------------------FUNCTIONS-----------------------

function [a, e, s_tilde] = my_lms(x, y, N, mu, normalized)
if nargin < 5
    normalized = false;
end
eps = 1e-05;
y = y(:);
x = x(:);

% initializing and zero padding
x       = [zeros(N-1, 1); x(:)];
L       = length(x);
a       = zeros(N, L-N+1);
e       = zeros(L-N+1, 1);
s_tilde = zeros(L-N+1, 1);


% main lms update loop
for n = 1:L-N+1
    xn = x(n:n+N-1);
    yn = y(n);

    % filter output and error calculation
    s_tilde(n) = a(:, n).'*xn;
    e(n)       = yn - s_tilde(n);

    % parameter update
    if normalized
        sigma     = sum(xn.^2)/N + eps;
        a(:, n+1) = a(:, n) + (mu/sigma)*e(n)*xn;
    else
        a(:, n+1) = a(:, n) + mu*e(n)*xn;
    end
end

a = flipud(a);
end
