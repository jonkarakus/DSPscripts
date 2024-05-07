%This MATLAB script, authored by Jon Karakus, demonstrates the implementation 
% of a feed-forward comb filter for noise cancellation. 
% It comprises two main sections: the first visualizes the amplitude response of
% the comb filter, while the second applies the comb filter and an LMS filter for noise cancellation.


%---------- Feed-forward Comb Filter--------------------- 

% %Define parameters
% M = 5; 
% b0 = 1; 
% g_values = [1, 0.75, 0.5, -1, -0.75, -0.5]; % Different values for g = bM / b0
% % Frequency vector for the plot
% omega = linspace(-pi, pi, 1000);
% figure;
% hold on;
% for g = g_values
%  bM = g * b0;
%  % Amplitude response |C(e^jω)|
%  amplitude_response = abs(b0 + bM * exp(-1j * omega * M));
% 
%  % Plot the amplitude response
%  plot(omega, amplitude_response);
% end
% title('Amplitude Response of Feed-forward Comb Filter');
% xlabel('\omega (radians)');
% ylabel('|C(e^{j\omega})|');
% legend(string(g_values), 'Location', 'Best');
% grid on;
% hold off;

%------------Applying Comb filter for Noise Cancellation------------
fs = 8000; 
M = round(fs / 60); 
b0 = 1; 
bM = -1; % Assuming bM is -1 to create a notch

[   d, Fs] = audioread('corrupted_speech.wav');

y = comb_filter(d, M, b0, bM);

D = fft(d);
Y = fft(y);
f = (0:length(D)/2-1)*(Fs/length(D));

figure;
subplot(2,1,1);
plot(f, 20*log10(abs(D(1:length(f)))));
title('Original Signal Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
subplot(2,1,2);
plot(f, 20*log10(abs(Y(1:length(f)))));
title('Filtered Signal Power Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%-------------Applying LMS filter for Noise Cancellation-----------
fs = 16000;
N = 16; 
mu = 0.01; 

[noise, fs_noise] = audioread('noise.wav');
[corrupted_speech, fs_speech] = audioread('corrupted_speech.wav');

if fs_noise ~= fs_speech
    noise = resample(noise, fs_speech, fs_noise);
    fs_noise = fs_speech;
end

reference_noise = noise / max(abs(noise));

input_signal = corrupted_speech / max(abs(corrupted_speech));

[a, error_signal, ~] = my_lms(reference_noise, input_signal, N, mu);


soundsc(error_signal, fs_speech);

t = (0:length(error_signal)-1) / fs_speech;
figure;
subplot(2, 1, 1);
plot(t, corrupted_speech);
title('Original Corrupted Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(t, error_signal);
title('Filtered Speech Signal (After LMS)');
xlabel('Time (s)');
ylabel('Amplitude');





%------------------FUNCTIONS--------------
function y = comb_filter(d, M, b0, bM)
    % d is the input signal, M is the delay length, b0 and bM are filter coefficients
    
    N = length(d);
    y = zeros(N, 1); % Initialize the output signal
    for n = M+1:N
        % Implement the difference equation of the comb filter
        y(n) = b0 * d(n) + bM * d(n-M);
    end
end

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



