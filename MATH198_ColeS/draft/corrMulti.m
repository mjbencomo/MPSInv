function X = corrMulti(kernel,Y)
%CORRMULTI Compute a discrete multi-channel cross-correlation
%   X = CORRMULTI(KERNEL,Y) computes the multi-channel cross-correlation
%
%       X(:,m) = sum_n corrDisc(KERNEL(:,n,m),Y(:,n)).
%
%   If Y has size K-by-N, then KERNEL must have size
%   (2*K-1)-by-N-by-M, where
%
%       KERNEL(:,n,m)
%
%   is the cross-correlation kernel mapping input channel m to output channel n.
%   The output X has size K-by-M.
%
%   For one output channel and multiple input channels, KERNEL must retain
%   its singleton output dimension; for example, its size should be
%   (2*K-1)-by-1-by-M rather than (2*K-1)-by-M.

arguments
    kernel {mustBeNumeric,mustBeNonempty,mustBeReal}
    Y      {mustBeNumeric,mustBeNonempty,mustBeReal}
end

if ~ismatrix(Y)
    error('corrMulti:InvalidInputDimensions', ...
        'Y must be a two-dimensional matrix with one column per output channel.');
end

if ndims(kernel) > 3
    error('corrMulti:InvalidKernelDimensions', ...
        'kernel must have size (2*K-1)-by-M-by-N.');
end

K = size(Y,1);
N = size(kernel,2);
M = size(kernel,3);

if size(kernel,1) ~= 2*K-1
    error('corrMulti:TimeSizeMismatch', ...
        'The first dimension of kernel must equal 2*size(X,1)-1.');
end

if size(Y,2) ~= N
    error('corrMulti:InputSizeMismatch', ...
        'The number of columns in Y must equal size(kernel,3).');
end

X = zeros(K,M,'like',Y);

for m = 1:M
    for n = 1:N
        X(:,m) = X(:,m) + corrDisc(kernel(:,n,m),Y(:,n));
    end
end
end
