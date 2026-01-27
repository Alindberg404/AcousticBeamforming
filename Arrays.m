clear all; clc;
%%
N = 5;
d = 0.05;
x = ((0:N-1) - (N-1)/2).' * d;
f = 1000;
c = 343; 
k = 2*pi*f/c;
theta0 = deg2rad(90);              %TARGETANGLE
phi_deg = 0:0.25:180;
phi = deg2rad(phi_deg);  
%%
SF = exp(-1j*k*x.*cos(theta0));    %SF - STEERFACTOR
TF = taper_uniform(N);             %TF - TAPERFACTOR
W = TF .* SF;                      %W - WEIGHTS

LV = exp(1j * k * (x * cos(phi))); %LV - LOOKING VECTOR

AF = W.' * LV;
AF_mag = abs(AF);
AF_mag = AF_mag / max(AF_mag);  % normalize

figure;
polarplot(phi, AF_mag, 'LineWidth', 2);
title('Array Factor |A(\theta)|');

%%
freqs = 200:20:4000;
B = zeros(numel(freqs), numel(phi));
TF = taper_uniform(N);

for fi = 1:numel(freqs)
    f = freqs(fi);
    k = 2*pi*f/c;
    SF = exp(-1j*k*x.*cos(theta0));
    W = TF .* SF;                      %W - WEIGHTS
    LV = exp(1j * k * (x * cos(phi))); %LV - LOOKING VECTOR

    AF = W.' * LV;
    AF_mag = abs(AF);
    B(fi,:) = AF_mag / max(AF_mag + eps);
end
figure;
imagesc(phi_deg, freqs/1000, B);   % x=phi in degrees, y=f in kHz
axis xy;
%
xlabel('\phi (deg)  (90° = broadside)');
ylabel('Frequency (kHz)');
colormap parula;
colorbar;
%%
a=taper_uniform(N)
N=10
b=taper_binomial(N)
%%

N_levels = [5, 7];
d_levels = [0.10, 0.02];

x = make_nested_array(N_levels, d_levels);

figure;
stem(x, ones(size(x)), 'filled');
xlabel('x [m]');
title('Nested linear array');
grid on;

%%
function x = make_nested_array(N_levels, d_levels)
    assert(numel(N_levels) == numel(d_levels), ...
        'N_levels and d_levels must have same length');

    x_all = [];

    for k = 1:numel(N_levels)
        N = N_levels(k);
        d = d_levels(k);

        xk = ((0:N-1) - (N-1)/2) * d;
        x_all  = [x_all, xk]; %#ok<AGROW>
    end
    x = unique(x_all);
    x = sort(x);
end

function w = taper_uniform(N)
w = ones(N,1);
w = w ./ sum(w);
end

function w = taper_binomial(N)
    k = 0:(N-1);
    w = arrayfun(@(kk) nchoosek(N-1, kk), k).';
    w = double(w);
    w = w ./ sum(w);
end