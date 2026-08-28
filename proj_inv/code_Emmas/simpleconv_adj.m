function w = simpleconv_adj(d, g)
% Computes cross correlation via FFT. Equivalent to convolution with
% conjugate g
%
% INPUTS:
% w - data vector (n x 1)
% g - Green's function (m x 1)
%
% OUTPUT:
% d - source vector (n x 1)

    n = length(d);
    m = length(g);
    L = n + m - 1;

    % Zero-padding
    g_ext = [g; zeros(L-m,1)];
    d_ext = [d; zeros(L-n,1)];

    % FFT convolution
    w_full = ifft(fft(d_ext).*conj(fft(g_ext)) );

    % Return first n entries
    w = real(w_full(1:n));
end

