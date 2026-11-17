function v = fractionalDerivative(u, dt, alpha, opts)
%FRACTIONALDERIVATIVE Apply the Grunwald-Letnikov derivative.
%   V = fractionalDerivative(U, DT, ALPHA) acts along the first nonsingleton
%   dimension of U. DT is the positive uniform time step; ALPHA >= 0.
%   ALPHA = 0 returns U. Inputs and outputs use double precision.
%
%   V = fractionalDerivative(..., Dimension=DIM) selects the time dimension.
%   V = fractionalDerivative(..., Transpose=true) applies the matrix transpose.
%   Real coefficients make this also the adjoint for complex-valued data.
%   Channels in all other dimensions are processed independently.
%
%   Samples preceding the first input sample are zero. The first output
%   is DT^(-ALPHA)*U(1), not an imposed zero initial value.
%   The sum at sample ell uses lags 0:ell-1 (MATLAB indexing).
%   Weights are generated recursively, without gamma functions.
%   No dense matrix is formed. Work is O(N*L) per channel, where L=N
%   except for an integer derivative, for which L=min(N,ALPHA+1).
%
%   Example:
%       t = (0:0.01:1).';
%       v = fractionalDerivative(sin(t), 0.01, 0.5);

arguments
    u double {mustBeFinite}
    dt (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
    alpha (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
    opts.Dimension double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive} = []
    opts.Transpose (1,1) logical = false
end

if ~isempty(opts.Dimension) && ~isscalar(opts.Dimension)
    error('fractionalDerivative:InvalidDimension', ...
        'Dimension must be empty or a positive integer scalar.');
end
v = u;
if isempty(u) || alpha == 0
    return
end

dim = opts.Dimension;
if isempty(dim)
    dim = find(size(u) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end
n = size(u, dim);
if alpha == fix(alpha)
    lastLag = min(n - 1, alpha);
else
    lastLag = n - 1;
end
weights = ones(lastLag + 1, 1);
for k = 1:lastLag
    weights(k + 1) = weights(k) * (k - 1 - alpha) / k;
end
weights = dt^(-alpha) * weights;

% For a lower triangular Toeplitz matrix A, A.' = J*A*J,
% where J reverses the sample order along the time dimension.
if opts.Transpose
    u = flip(u, dim);
end
v = filter(weights, 1, u, [], dim);
if opts.Transpose
    v = flip(v, dim);
end
end
