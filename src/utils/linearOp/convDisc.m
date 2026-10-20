function y = convDisc(kernel,x)
%CONVDISC Discrete convolution via FFT-based circular convolution.
%   y = CONVDISC(kernel,x) computes the discrete convolution of the real,
%   numeric vectors kernel and x by converting it to a circular
%   convolution problem and evaluating that with FFTs.
%
%   Assumptions:
%   - kernel and x must be nonempty real vectors
%   - numel(kernel) must equal 2*numel(x)-1
%   - the output y has the same orientation as x
%
%   The implementation forms zero-padded, shifted vectors and calls
%   CONVCIRC, which computes the circular convolution using FFTs.

arguments
    kernel {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
    x      {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeReal}
end

n = numel(x);

if numel(kernel) ~= 2*n - 1
    error('convDisc:SizeMismatch', ...
        'The kernel must have length 2*numel(x)-1.');
end

% Work with column vectors internally.
outputIsRow = isrow(x);
kernel = kernel(:);
x = x(:);

% Step 1: construct the modified vectors.
kernelZ = [kernel(n:end); ...
    zeros(1,1,'like',kernel); ...
    kernel(1:n-1)];

xZ = [x; zeros(n,1,'like',x)];

% Step 2: compute their circular convolution
yZ = convCirc(kernelZ,xZ);

% Step 3: extract the first n entries.
y = yZ(1:n);

% Restore the orientation of y.
if outputIsRow
    y = y.';
end
end