function d = simpleconv(w, g)
% Computes convolution using FFT
%
% INPUTS:
% w - source vector (n x 1)
% g - Green's function kernel (m x 1)
%
% OUTPUT:
% d - truncated linear convolution (n x 1)
    n = length(w);
    m = length(g);
    L = n + m - 1;

    % Zero-padding
    g_ext = [g; zeros(L-m,1)];
    w_ext = [w; zeros(L-n,1)];

    % FFT convolution
    d_full = ifft( fft(g_ext) .* fft(w_ext) );

    % Return all entries
    d = real(d_full(1:n));
end

