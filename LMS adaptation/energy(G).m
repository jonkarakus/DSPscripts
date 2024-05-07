%Code authored by Jon Karakus to investigate energy and standard deviation
%for Gaussian noise and speech files. 

[speech, fs] = audioread('speech.wav');

% Compute the energy G for the entire file
G = sum(abs(speech).^2);

% Get the duration of the signal
duration = length(speech) / fs;

% Calculate the number of samples
num_samples = round(duration * fs);

% Generate zero-mean Gaussian noise with the same parameters
std_dev = sqrt(G / num_samples); % Standard deviation for Gaussian noise
gaussian_noise = normrnd (0, std_dev, num_samples, 1);

% Write the Gaussian noise to gauss.wav file
audiowrite('gauss.wav', gaussian_noise, fs);

% Verify that gauss.wav has the same energy as speech.wav
[gauss, ~] = audioread('gauss.wav');
G_gauss = sum(abs(gauss).^2);

% Display the value of G and the standard deviation used
disp(['Energy G for speech.wav: ', num2str(G)]);
disp(['Standard deviation for Gaussian noise: ', num2str(std_dev)]);
disp(['Energy G for gauss.wav: ', num2str(G_gauss)]);