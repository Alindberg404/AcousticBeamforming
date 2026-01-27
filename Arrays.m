clear; clc;

N = 5;
d = 0.05;
x = ((0:N-1) - (N-1)/2) * d;   % element positions (row)
f = 1000;
c = 343;
k = 2*pi*f/c;

phi_deg = 0:0.25:180;
phi = deg2rad(phi_deg);

w = taper_uniform(N);          % make sure this returns a ROW vector 1xN

% Steering vector for every phi -> N x numAngles
SV = exp(1j*k*(x.').*cos(phi));    % x' is Nx1, cos(phi) is 1xM => NxM

AF = w' * SV;                   % (1xN)*(NxM) = 1xM
mag = abs(AF);
mag = mag / max(mag);          % normalize (optional)

figure;
polarplot(phi, mag, "LineWidth", 2);


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