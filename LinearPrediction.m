% In this MATLAB script authored by Jon Karakus, an AutoRegressive Moving Average (ARMA) 
% signal is analyzed and compared through two key processes. Initially, the script computes 
% the mean-square prediction error across filter orders to assess the performance. Subsequently, 
% it synthesizes a new signal and juxtaposes its autocorrelation 
% with that of the original, providing a view of signal fidelity and
% synthesis. 


%--ARMA signal ploted as mean-square prediction error against order------

x = [1, 0.2, 0.3]; %AR filter coefficients in denominator
b = [1, -0.6]; %AR filter coefficients in numerator
arma = filter(b,x, randn(8000,1)); %pass AR coefficients and zero mean and unit variance 
gaussian white noise into LTI filter
mse = zeros(10,1); %initialize mean-square error array
for p = 1:10 %loop over p
 [a,g] = lpc(arma,p); %compute prediction error filter coefficients and mean-square 
prediction error for ARMA signal
 mse(p) = g;
 fprintf('order: %d: mean square prediction error: %f\n',p,mse(p));
end
figure; %plot mean-square prediction error for filter orders from 1 to 10
plot(1:10,mse);
xlabel('filter order (p)');
ylabel('mean-square prediction error');
title('mean-square prediction error vs. filter order');
grid on;


%--------------AutoCorrelated Synthesized Fucntions------------------

x = [1, 0.2, 0.3]; %AR filter coefficients in denominator
b = [1, -0.6]; %AR filter coefficients in numerator
noise = randn(8000,1);
arma = filter(b,x,noise); %input white gaussian noise with 0 mean and unit variance into filter
p = 7; %set filter order to 7
[a,g] = lpc(arma,p); %determine prediction error filter coefficients and mean-square prediction 
%error (used for value of DC gain constant G in the causal all-pole system used to synthesize the signal)
output = filter(g,a,noise); %inputting white noise into all-pole filter to synthesize the signal
autocorr1 = xcorr(arma); %autocorrelation array for the original signal (from question 2)
autocorr2 = xcorr(output); %autocorrelation array for the synthesized signal
figure; %plot autocorrelation values for the original and synthesized signals on same graph for comparison
plot(autocorr1, 'b', 'LineWidth', 3);
hold on;
plot(autocorr2, 'r', 'LineWidth', 1);
hold off;
xlabel('Lag');
ylabel('Autocorrelation');
title('Autocorrelation Comparison');
legend('original signal', 'modified signal');
grid on;
numMatchingTerms = sum(abs(autocorr1 - autocorr2) < 1e-7); %find how many autocorrelation terms 
match
fprintf('Number of matching autocorrelation terms: %d\n', numMatchingTerms);



