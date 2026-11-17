classdef MPSPreconditionerOp < InvLinearOp
    %MPSPRECONDITIONEROP Physics-based preconditioner B = Q'*Q, Sec. 9.6.
    %   B = MPSPreconditionerOp(c,d,s,Nt,dt)
    %
    %   c   Positive constant wave speed.
    %   d   Spatial dimension: 1, 2, or 3.
    %   s   M-by-d array; row m is the nonnegative integer multi-index s_m.
    %   Nt  Positive number of time samples per channel.
    %   dt  Positive uniform time step.
    %
    %   For r_m = sum(s(m,:)) and alpha_m = (d-1)/2 + r_m,
    %       Q_m = c^(-r_m) * D^(alpha_m).
    %   fwd applies Q'*Q; inv applies Q^(-1)*Q^(-T).
    %   adj = fwd and invAdj = inv since B is real symmetric positive
    %   definite under the fractional functions' zero-prehistory convention.
    %
    %   All actions take a column vector of length Nt*M. Channel samples
    %   are stacked using x = X(:), where X is Nt-by-M.
    %   Uses Euclidean adjoints and the existing fractionalDerivative and
    %   fractionalIntegral functions, including their Transpose option.
    %
    %   Expected superclass constructor: InvLinearOp(n).
    %   If your superclass takes two dimensions, change the superclass
    %   constructor call below to obj@InvLinearOp(n,n).

    properties (SetAccess = private)
        c (1,1) double
        d (1,1) double
        s (:,:) double
        Nt (1,1) double
        dt (1,1) double
        M (1,1) double
        beta (1,1) double
        orders (:,1) double
        alpha (:,1) double
    end

    methods
        function obj = MPSPreconditionerOp(c,d,s,Nt,dt)
            arguments
                c (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
                d (1,1) double {mustBeMember(d,[1,2,3])}
                s (:,:) double ...
                    {mustBeReal,mustBeFinite,mustBeInteger,mustBeNonnegative}
                Nt (1,1) double ...
                    {mustBeReal,mustBeFinite,mustBeInteger,mustBePositive}
                dt (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
            end
            if isempty(s) || size(s,2) ~= d
                error('MPSPreconditionerOp:InvalidMultiIndices', ...
                    's must be a nonempty M-by-d array of multi-indices.');
            end
            n = Nt*size(s,1);
            obj@InvLinearOp(n);
            obj.c = c;
            obj.d = d;
            obj.s = s;
            obj.Nt = Nt;
            obj.dt = dt;
            obj.M = size(s,1);
            obj.beta = (d-1)/2;
            obj.orders = sum(s,2);
            obj.alpha = obj.beta + obj.orders;
        end

        function y = fwd(obj,x)
            % Apply B_m = c^(-2*r_m) * D_m' * D_m.
            obj.validateDomainVector(x);
            X = reshape(double(x),obj.Nt,obj.M);
            Y = zeros(size(X),'like',X);
            for m = 1:obj.M
                z = fractionalDerivative(X(:,m),obj.dt,obj.alpha(m), ...
                    Dimension=1);
                z = fractionalDerivative(z,obj.dt,obj.alpha(m), ...
                    Dimension=1,Transpose=true);
                Y(:,m) = obj.c^(-2*obj.orders(m))*z;
            end
            y = Y(:);
        end

        function x = adj(obj,y)
            obj.validateRangeVector(y);
            x = obj.fwd(y);
        end

        function x = inv(obj,b)
            % Apply B_m^(-1) = c^(2*r_m) * I_m * I_m'.
            % Apply the transposed integral FIRST; order matters.
            obj.validateRangeVector(b);
            R = reshape(double(b),obj.Nt,obj.M);
            X = zeros(size(R),'like',R);
            for m = 1:obj.M
                z = fractionalIntegral(R(:,m),obj.dt,obj.alpha(m), ...
                    Dimension=1,Transpose=true);
                z = fractionalIntegral(z,obj.dt,obj.alpha(m), ...
                    Dimension=1);
                X(:,m) = obj.c^(2*obj.orders(m))*z;
            end
            x = X(:);
        end

        function x = invAdj(obj,b)
            obj.validateDomainVector(b);
            x = obj.inv(b);
        end

        function B = matrixRep(obj)
            % Explicit dense representation for small problems only.
            % Forward and inverse actions do not call this method.
            blocks = cell(1,obj.M);
            for m = 1:obj.M
                D = fractionalDerivative(eye(obj.Nt),obj.dt,obj.alpha(m), ...
                    Dimension=1);
                blocks{m} = obj.c^(-2*obj.orders(m))*(D'*D);
            end
            B = blkdiag(blocks{:});
        end
    end
end
