function x = corrDisc(kernel,y)
%CORRDISC Discrete cross-correlation via FFT-based circular cross-correlation.
%   x = CORRDISC(kernel,y) computes the discrete cross-correlation of the real,
%   numeric vectors kernel and y by converting it to a circular
%   cross-correlation problem and evaluating that with FFTs.
%
%   Assumptions:
%   - kernel and y must be nonempty real vectors
%   - numel(kernel) must equal 2*numel(y)-1
%   - the output x has the same orientation as y
%
%   The implementation forms zero-padded, shifted vectors and calls
%   CORRCIRC, which computes the circular cross-correlation using FFTs.

arguments
    kernel {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
    y      {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
end

n = numel(y);

if numel(kernel) ~= 2*n - 1
    error('corrDisc:SizeMismatch', ...
        'The kernel must have length 2*numel(y)-1.');
end

% Work with column vectors internally.
outputIsRow = isrow(y);
kernel = kernel(:);
y = y(:);

% Step 1: construct the modified vectors.
kernelZ = [kernel(n:end); ...
    zeros(1,1,'like',kernel); ...
    kernel(1:n-1)];

yZ = [y; zeros(n,1,'like',y)];

% Step 2: compute their circular convolution
xZ = corrCirc(kernelZ,yZ);

% Step 3: extract the first n entries.
x = xZ(1:n);

% Restore the orientation of x.
if outputIsRow
    x = x.';
end
end