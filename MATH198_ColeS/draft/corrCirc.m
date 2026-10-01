function x = corrCirc(kernel,y)
%CORRCIRC Circular cross-correlation of two real vectors using FFTs.
%   x = CORRCIRC(kernel,y) returns the circular cross-correlation of the real,
%   numeric vectors kernel and y. Both inputs must be nonempty vectors of
%   the same size.
%
%   The result is computed with FFT-based multiplication in the frequency
%   domain.

arguments
    kernel {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
    y      {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
end

if ~isequal(size(kernel),size(y))
    error('corrCirc:SizeMismatch', ...
        'kernel and y must have the same size.');

end

x = ifft(conj(fft(kernel)) .* fft(y), 'symmetric');
end