clear; clc; close all;
load handel.mat


% freqs : column vector [Hz] (e.g. (f1:fr:f2).')
% H     : column vector complex same length as freqs, where H(fi)=AF at your chosen angle

freqs = (300:20:3400).';
H = ones(size(freqs));  % <-- replace with your extracted complex H(fi)

N = length(y);
Y = fft(y);

% FFT frequency bins for the positive half (0..Fs/2)
K = floor(N/2) + 1;
f_half = (0:K-1).' * (Fs/N);

% Interpolate complex H onto these bins
% Outside your defined band -> 0 (so you only apply what you know)
H_half = interp1(freqs, H, f_half, 'linear', 0);

% Optional: normalize overall gain so it doesn't clip / get too quiet
H_half = H_half ./ (max(abs(H_half)) + eps);

% Build full conjugate-symmetric response H_full (length N)
H_full = zeros(N,1);
H_full(1:K) = H_half;

if mod(N,2)==0
    % even length: bins are [0 .. Fs/2] then negative freqs
    H_full(K+1:end) = conj(H_half(K-1:-1:2));
else
    % odd length: no exact Nyquist bin
    H_full(K+1:end) = conj(H_half(K:-1:2));
end

% Apply and transform back
Y = Y .* H_full;
x_filt = ifft(Y, 'symmetric');

% Normalize for listening
x_filt = x_filt ./ (max(abs(x_filt)) + eps);

%% --- Listen ---
soundsc(y, Fs);       % original
pause(3);
soundsc(x_filt, Fs);  % "heard through" your array response at that angle

%% --- Optional: plot magnitude response you applied ---
figure;
plot(f_half, 20*log10(abs(H_half)+eps), 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('|H| (dB)');
title('Applied angle response (interpolated to FFT bins)');
