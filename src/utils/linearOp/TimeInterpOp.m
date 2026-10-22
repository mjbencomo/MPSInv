classdef TimeInterpOp < LinearOp
    % Piecewise-linear time interpolation, independently on each channel.
    %
    % P = TimeInterpOp(tInput,tOutput,NumChannels=C) maps X(:) to Y(:),
    % where X is numel(tInput)-by-C and Y is numel(tOutput)-by-C.
    % tInput must increase strictly; tOutput may be in any order, including
    % repeated times. Both inputs are nonempty, finite, real time vectors.
    % Grid1D objects can be used by passing their .pts vectors.
    %
    % OutsideDomain="error" (default) rejects times outside the input span.
    % OutsideDomain="zero" gives zero at those times. A singleton input
    % grid is supported only at its exact time, with this same outside rule.
    %
    % adj is the Euclidean/Hermitian adjoint, not reverse interpolation.
    % There is no time-step factor in interpolation or its adjoint.
    % fwd and adj return full vectors; weights and matrixRep remain sparse.

    properties (SetAccess = private)
        tInput (:,1) double
        tOutput (:,1) double
        numChannels (1,1) double
        outsideDomain (1,1) string
        weights (:,:) double  % Sparse single-channel interpolation matrix.
    end

    methods
        function obj = TimeInterpOp(tInput,tOutput,options)
            arguments
                tInput double {mustBeReal,mustBeFinite,mustBeNonempty}
                tOutput double {mustBeReal,mustBeFinite,mustBeNonempty}
                options.NumChannels (1,1) double ...
                    {mustBeInteger,mustBePositive,mustBeFinite} = 1
                options.OutsideDomain (1,1) string ...
                    {mustBeMember(options.OutsideDomain,["error","zero"])} = "error"
            end
            if ~isvector(tInput) || ~isvector(tOutput)
                error('TimeInterpOp:InvalidTimeVector', ...
                    'Input and output times must be vectors.');
            end
            tInput = tInput(:);
            tOutput = tOutput(:);
            if any(diff(tInput) <= 0) || any(~isfinite(diff(tInput)))
                error('TimeInterpOp:InvalidInputTimes', ...
                    'Input times must increase strictly with finite intervals.');
            end

            K = numel(tInput);
            L = numel(tOutput);
            outside = tOutput < tInput(1) | tOutput > tInput(end);
            if any(outside) && options.OutsideDomain == "error"
                error('TimeInterpOp:OutsideDomain', ...
                    'Output times must lie within the input time interval.');
            end

            rows = find(~outside);
            if K == 1
                W = sparse(rows,ones(size(rows)),ones(size(rows)),L,1);
            else
                % discretize assigns the final endpoint to the final bin.
                left = discretize(tOutput(rows),tInput);
                theta = (tOutput(rows)-tInput(left))./ ...
                    (tInput(left+1)-tInput(left));
                W = sparse([rows;rows],[left;left+1], ...
                    [1-theta;theta],L,K);
            end

            obj@LinearOp(K*options.NumChannels,L*options.NumChannels);
            obj.tInput = tInput;
            obj.tOutput = tOutput;
            obj.numChannels = options.NumChannels;
            obj.outsideDomain = options.OutsideDomain;
            obj.weights = W;
        end

        function y = fwd(obj,x)
            obj.validateDomainVector(x);
            X = reshape(x,numel(obj.tInput),obj.numChannels);
            Y = obj.weights*X;
            y = full(Y(:));
        end

        function x = adj(obj,y)
            obj.validateRangeVector(y);
            Y = reshape(y,numel(obj.tOutput),obj.numChannels);
            X = obj.weights'*Y;
            x = full(X(:));
        end

        function A = matrixRep(obj)
            A = kron(speye(obj.numChannels),obj.weights);
        end
    end
end
