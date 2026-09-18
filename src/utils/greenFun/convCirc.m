function y = convCirc(kernel,x)
%CONVCIRC Circular convolution of two real vectors using FFTs.
%   y = CONVCIRC(kernel,x) returns the circular convolution of the real,
%   numeric vectors kernel and x. Both inputs must be nonempty vectors of
%   the same size.
%
%   The result is computed with FFT-based multiplication in the frequency
%   domain.

arguments
    kernel {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
    x      {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
end

if ~isequal(size(kernel),size(x))
    error('convCirc:SizeMismatch', ...
        'kernel and x must have the same size.');

end

y = ifft(fft(kernel) .* fft(x), 'symmetric');
end