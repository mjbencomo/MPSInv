function d = mcconv(w, g)
% Multichannel convolution with one output channel (need to adjust to take
% multiple output channels)
%
% INPUTS:
% w  - (n x C_in)
% G  - (m x C_in)
%
% OUTPUT:
% d  - (n x 1)

    [n, C] = size(w);
    [m, ~] = size(g);

    L = n + m - 1;
    d_full = zeros(L,1);

    for c = 1:C
        wc = w(:,c);
        gc = g(:,c);

        d_full = d_full + ifft(fft(wc, L) .* fft(gc, L));
    end

    d = real(d_full(1:n));
end