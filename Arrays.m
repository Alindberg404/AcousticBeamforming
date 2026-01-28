clear all; clc;close all;
%%
close all;
N=5;d=0.05;theta0=90-20;f=1000;f1=300;fr=20;f2=3400;

%plot_arrayfactor_polar(x_array, f, theta0, TF)

plot_arrayfactor_heatmap(N,d,'simple','uni',f1,fr,f2,theta0,'fix')
plot_arrayfactor_heatmap(N,d,'simple','uni',f1,fr,f2,theta0,'fdep')


%%
x_array = make_nested_array([5,8],[0.5,0.02]);
visualize_nested_array(x_array)

%%
%#####################################################################
%#####################################################################
function plot_arrayfactor_heatmap(N, d, array_mode, taper_mode, f1,fr,f2, theta0_deg, w_mode)
    c=343;
    freqs = f1:fr:f2;
    theta0 = deg2rad(theta0_deg);
    phi_deg = 0:0.25:180;
    phi = deg2rad(phi_deg); 
    
    B = zeros(numel(freqs), numel(phi));

    if strcmp(array_mode, 'simple')
        x_array=make_simple_array(N,d);
    end
    if strcmp(array_mode, 'nested')
        %x_array=make_nested_array(N,d);
    end
    
    if strcmp(taper_mode, 'uni')
        TF=taper_uniform(N);
    end

    if strcmp(taper_mode, 'bin')
        TF=taper_binomial(N);
    end
    

    if strcmp(w_mode, 'fix')
        f_ref = mean(freqs);
        k_ref = 2*pi*f_ref/c;
        SF = exp(-1j*k_ref*x_array.*cos(theta0));
        W = TF .* SF; 
        for fi = 1:numel(freqs)
            f = freqs(fi);
            k = 2*pi*f/c;
            LV = exp(1j * k * (x_array * cos(phi))); %LV - LOOKING VECTOR
            AF = W.' * LV;
            AF_mag = abs(AF);
            B(fi,:) = AF_mag / max(AF_mag + eps);
        end
    end
    
    if strcmp(w_mode, 'fdep')
        for fi = 1:numel(freqs)
            f = freqs(fi);
            k = 2*pi*f/c;
            SF = exp(-1j*k*x_array.*cos(theta0));
            W = TF .* SF;                            %W - WEIGHTS
            LV = exp(1j * k * (x_array * cos(phi))); %LV - LOOKING VECTOR
        
            AF = W.' * LV;
            AF_mag = abs(AF);
            B(fi,:) = AF_mag / max(AF_mag + eps);
        end
    end

    figure;
    imagesc(phi_deg, freqs/1000, B);   % x=phi in degrees, y=f in kHz
    axis xy;
    title(sprintf(['Array Factor |A(\\phi,f)| — %s array, %s taper, %s weights\n' ...
                   'N = %d, d = %.2f m, \\theta_0 = %d°'], ...
                   array_mode, taper_mode, w_mode, ...
                   N, d, theta0_deg));    xlabel('\phi (deg)  (90° = broadside)');
    ylabel('Frequency (kHz)');
    colormap parula;
    colorbar;
    hold on;
    xline(theta0_deg, 'b--', 'LineWidth', 1.0);
    hold off;
    
    plotdb=0;
    if plotdb == 1
        figure;
        lim=-20;
        BdB = 20*log10(B + eps);
        BdB(BdB < lim) = lim;      % clip dynamic range
        imagesc(phi_deg, freqs/1000, BdB); axis xy;
        clim([lim 0]); 
        colormap parula;
        colorbar;
        title(sprintf('AF dB mode: %s, \\theta_0 = %d°', w_mode, theta0_deg));
        hold on;
        xline(theta0_deg, 'b--', 'LineWidth', 1.0);
        hold off;
    end
end

function plot_arrayfactor_polar(x_array, f, theta0_deg, TF)
    theta0 = deg2rad(theta0_deg);
    phi_deg = 0:0.25:180;
    phi = deg2rad(phi_deg); 
    
    c = 343; 
    k = 2*pi*f/c;
    SF = exp(-1j*k*x_array.*cos(theta0));    %SF - STEERFACTOR
    W = TF .* SF;                      %W - WEIGHTS
    
    LV = exp(1j * k * (x_array * cos(phi))); %LV - LOOKING VECTOR
    
    AF = W.' * LV;
    AF_mag = abs(AF);
    AF_mag = AF_mag / max(AF_mag);  % normalize
    
    figure;
    polarplot(phi, AF_mag, 'LineWidth', 2);
    title('Array Factor |A(\theta)|');
end

function x = make_simple_array(N,d)
    x = ((0:N-1) - (N-1)/2).' * d;
end

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

function visualize_nested_array(x)
    figure;
    stem(x, ones(size(x)), 'filled');
    xlabel('x [m]');
    title('Nested linear array');
    grid on;
end