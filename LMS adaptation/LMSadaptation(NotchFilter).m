%This MATLAB script, authored by Jon Karakus, demonstrates adaptive noise cancellation (ANC) employing the 
%Least Mean Squares (LMS) algorithm. It aims to attenuate tonal noise from a primary signal containing both 
%noise and a content-carrying signal. The script iterates over different content signal frequencies, applies 
%the LMS algorithm, and calculates the relative power of the output. Additionally, it analyzes the notch width 
%around the reference frequency to assess noise suppression effectiveness.


%---------------LMS adaptation with tonal input(notch)-----------------
fs = 16000; 
duration = 5; 
amp = 0.3; 
f_noise = 1000; 
N = 16;
mu = 0.159; %learning rate

% Generate time vector
t = 0:1/fs:duration-1/fs;

% Generate the noise source signal
noise_source = amp * sin(2 * pi * f_noise * t);

frequencies = 500:1500; 
relative_powers = zeros(length(frequencies), 1);

% LMS adaptation loop
for i = 1:length(frequencies)
    f_content = frequencies(i);

    content_signal = amp * sin(2 * pi * f_content * t);

    primary_signal = content_signal + noise_source;
    % Reference signal (same as noise source, as H[z] = J[z] = 0)
    reference_signal = noise_source;

    % Apply the LMS algorithm
    [a, e, ~] = my_lms(reference_signal, primary_signal, N, mu);

    % Calculate the relative power of the output
    relative_powers(i) = sum(e.^2) / sum(primary_signal.^2);
end

figure;
plot(frequencies, 10*log10(relative_powers));
title('Relative Power of Output vs. Content-Carrying Signal Frequency');
xlabel('Frequency of Content-Carrying Signal (Hz)');
ylabel('Relative Power (dB)');
grid on;


%----------------------------Observation---------------
% % Load the .wav file
% [y, Fs] = audioread('corrupted_speech.wav');
% 
% % Duration of the audio file in seconds
% duration = length(y) / Fs;
% 
% %FFT 
% L = length(y);
% Y = fft(y);
% f = Fs*(0:(L/2))/L;
% 
% % Plot the waveform
% plot(f, abs(Y(1:L/2+1)));
% title('Frequency of Speech Audio');
% xlabel('Frequency');
% ylabel('Amplitude');

%--------------Notch Width around reference frequency--------------------

L = length(e); 
E = fft(e)
P2 = abs(E/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

notchPoint = max(10*log10(P1));
notchIndex = find(10*log10(P1) == notchPoint);
P1_3dB = notchPoint - 3;
above3dBIndices = find(10*log10(P1) > P1_3dB)

aboveLeft = above3dBIndices(above3dBIndices < notchIndex);
aboveRight = above3dBIndices(above3dBIndices > notchIndex);
if isempty(aboveLeft) || isempty(aboveRight)
    disp('Cannot determine notch width at -3 dB threshold.');
else
    f_low = f(max(aboveLeft));
    f_high = f(min(aboveRight));
    notchWidth_3dB = f_high - f_low;

    % Display the notch width
    disp(['Notch Width at -3 dB: ', num2str(notchWidth_3dB), ' Hz']);
    
    % Optional: Plot for visual verification
    figure;
    plot(f, 10*log10(P1)); hold on;
    plot([f_low f_high], [P1_3dB P1_3dB], 'r--'); % -3 dB line
    title('Error Signal Frequency Response and Notch Width at -3 dB');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    grid on; hold off;
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
