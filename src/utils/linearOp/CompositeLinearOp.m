classdef CompositeLinearOp < LinearOp
    % Composition in mathematical product order.
    %
    % C = CompositeLinearOp({A,B,D}) represents A o B o D:
    %   C.fwd(x)      = A.fwd(B.fwd(D.fwd(x)))
    %   C.adj(y)      = D.adj(B.adj(A.adj(y)))
    %   C.matrixRep() = A.matrixRep()*B.matrixRep()*D.matrixRep()
    %
    % Accepts a nonempty cell vector of scalar LinearOp objects, or a
    % nonempty LinearOp object vector. Use cells for mixed subclasses.
    % The rightmost operator acts first. Adjacent operators must satisfy
    % operators{k}.dimDomain == operators{k+1}.dimRange.
    %
    % Compatibility checks dimensions only. Callers must also ensure that
    % matching coordinates use the same grids, units, and channel ordering.
    % Forward/adjoint actions do not assemble matrix representations.

    properties (SetAccess = private)
        operators (1,:) cell
        numOperators (1,1) double
    end

    methods
        function obj = CompositeLinearOp(operators)
            narginchk(1,1);
            if isempty(operators) || ~isvector(operators)
                error('CompositeLinearOp:InvalidOperators', ...
                    'Provide a nonempty cell vector or LinearOp object vector.');
            end

            if iscell(operators)
                ops = reshape(operators,1,[]);
            elseif isa(operators,'LinearOp')
                ops = cell(1,numel(operators));
                for k = 1:numel(operators)
                    ops{k} = operators(k);
                end
            else
                error('CompositeLinearOp:InvalidOperators', ...
                    'Provide a cell vector or LinearOp object vector.');
            end

            for k = 1:numel(ops)
                if ~isa(ops{k},'LinearOp') || ~isscalar(ops{k})
                    error('CompositeLinearOp:InvalidOperator', ...
                        'Entry %d must be a scalar LinearOp object.',k);
                end
            end
            for k = 1:numel(ops)-1
                if ops{k}.dimDomain ~= ops{k+1}.dimRange
                    error('CompositeLinearOp:IncompatibleDimensions', ...
                        ['Operator %d (%s) has domain dimension %d, but ' ...
                         'operator %d (%s) has range dimension %d.'], ...
                        k,class(ops{k}),ops{k}.dimDomain, ...
                        k+1,class(ops{k+1}),ops{k+1}.dimRange);
                end
            end

            obj@LinearOp(ops{end}.dimDomain,ops{1}.dimRange);
            obj.operators = ops;
            obj.numOperators = numel(ops);
        end

        function y = fwd(obj,x)
            obj.validateDomainVector(x);
            y = x;
            for k = obj.numOperators:-1:1
                y = obj.operators{k}.fwd(y);
            end
        end

        function x = adj(obj,y)
            obj.validateRangeVector(y);
            x = y;
            for k = 1:obj.numOperators
                x = obj.operators{k}.adj(x);
            end
        end

        function A = matrixRep(obj)
            % Multiply component matrices, without sampling forward actions.
            A = obj.operators{end}.matrixRep();
            for k = obj.numOperators-1:-1:1
                A = obj.operators{k}.matrixRep()*A;
            end
        end
    end
end
