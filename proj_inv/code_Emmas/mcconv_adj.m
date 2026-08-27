function w = mcconv_adj(d, g)
% Multichannel cross-correlation (adjoint of convolution)
%
% INPUTS:
% d  - (n x 1)
% g  - (m x C_in)
%
% OUTPUT:
% w - (n x C_in)

    n = length(d);
    [m, C] = size(g);

    L = n + m - 1;

    w = zeros(n, C);

    for c = 1:C
        gc = g(:,c);
        tmp = ifft(fft(d,L) .* conj(fft(gc, L)));
        w(:,c) = real(tmp(1:n));
    end
end