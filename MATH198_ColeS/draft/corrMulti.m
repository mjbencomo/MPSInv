function Y = corrMulti(kernel,X)
%CORRMULTI Compute a discrete multi-channel cross-correlation
%   Y = CORRMULTI(KERNEL,X) computes the multi-channel cross-correlation
%
%       Y(:,n) = sum_m corrDisc(KERNEL(:,n,m),X(:,m)).
%
%   If X has size K-by-M, then KERNEL must have size
%   (2*K-1)-by-N-by-M, where
%
%       KERNEL(:,n,m)
%
%   is the cross-correlation kernel mapping input channel m to output channel n.
%   The output Y has size K-by-N.
%
%   For one output channel and multiple input channels, KERNEL must retain
%   its singleton output dimension; for example, its size should be
%   (2*K-1)-by-1-by-M rather than (2*K-1)-by-M.

arguments
    kernel {mustBeNumeric,mustBeNonempty,mustBeReal}
    X      {mustBeNumeric,mustBeNonempty,mustBeReal}
end

if ~ismatrix(X)
    error('corrMulti:InvalidInputDimensions', ...
        'X must be a two-dimensional matrix with one column per input channel.');
end

if ndims(kernel) > 3
    error('corrMulti:InvalidKernelDimensions', ...
        'kernel must have size (2*K-1)-by-N-by-M.');
end

K = size(X,1);
N = size(kernel,2);
M = size(kernel,3);

if size(kernel,1) ~= 2*K-1
    error('corrMulti:TimeSizeMismatch', ...
        'The first dimension of kernel must equal 2*size(X,1)-1.');
end

if size(X,2) ~= M
    error('corrMulti:InputSizeMismatch', ...
        'The number of columns in X must equal size(kernel,3).');
end

Y = zeros(K,N,'like',X);

for n = 1:N
    for m = 1:M
        Y(:,n) = Y(:,n) + corrDisc(kernel(:,n,m),X(:,m));
    end
end
end
