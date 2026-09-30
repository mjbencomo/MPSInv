classdef ConvMultiOp < LinearOp
    % ConvMultiOp implements discrete multi-channel convolution.
    %
    % If kernels has size (2*K-1)-by-N-by-M, then the operator maps M
    % input channels of length K to N output channels of length K. Domain
    % and range vectors use channel-block ordering:
    %
    %   x = X(:),  size(X) = [K,M]
    %   y = Y(:),  size(Y) = [K,N].

    properties (SetAccess = private)
        kernels
        scale (1,1) double = 1
        N_time (1,1) double
        N_input (1,1) double
        N_output (1,1) double
    end

    methods
        function obj = ConvMultiOp(kernels,options)
            arguments
                kernels {mustBeNumeric,mustBeNonempty,mustBeReal}
                options.Scale (1,1) double ...
                    {mustBeFinite,mustBeReal} = 1
            end

            if ndims(kernels) > 3
                error('ConvMultiOp:InvalidKernelDimensions', ...
                    'kernels must have size (2*K-1)-by-N-by-M.');
            end

            numberOfLags = size(kernels,1);
            if mod(numberOfLags,2) ~= 1
                error('ConvMultiOp:InvalidKernelTimeSize', ...
                    ['The first kernel dimension must be odd and equal ' ...
                     'to 2*K-1 for some positive integer K.']);
            end

            K = (numberOfLags+1)/2;
            N = size(kernels,2);
            M = size(kernels,3);

            obj@LinearOp(K*M,K*N);

            obj.kernels = kernels;
            obj.scale = options.Scale;
            obj.N_time = K;
            obj.N_input = M;
            obj.N_output = N;
        end

        function y = fwd(obj,x)
            % Apply scaled multi-channel convolution.
            obj.validateDomainVector(x);
            if ~isreal(x)
                error('ConvMultiOp:ComplexInput', ...
                    'The current convMulti implementation requires real data.');
            end

            X = reshape(x,obj.N_time,obj.N_input);
            Y = obj.scale*convMulti(obj.kernels,X);
            y = Y(:);
        end

        function x = adj(obj,y)
            % Apply the Euclidean adjoint using multi-channel correlation.
            obj.validateRangeVector(y);
            if ~isreal(y)
                error('ConvMultiOp:ComplexInput', ...
                    'The current convMulti implementation requires real data.');
            end

            Y = reshape(y,obj.N_time,obj.N_output);
            X = obj.scale*corrMulti(adjointKernels,Y); %MB: Need to implement cross correlation
            x = X(:);
        end

        function A = matrixRep(obj)
            % Assemble the block-Toeplitz matrix representation.
            K = obj.N_time;
            N = obj.N_output;
            M = obj.N_input;

            A = zeros(obj.dimRange,obj.dimDomain,'like',obj.kernels);

            for n = 1:N
                rows = (n-1)*K+(1:K);

                for m = 1:M
                    columns = (m-1)*K+(1:K);
                    kernel = obj.kernels(:,n,m);

                    firstColumn = kernel(K:end);
                    firstRow = kernel(K:-1:1).';

                    A(rows,columns) = obj.scale*toeplitz( ...
                        firstColumn,firstRow);
                end
            end
        end
    end
end
