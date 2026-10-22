classdef AjointOfLinearOp < LinearOp
    % Operator wrapper representing the Euclidean/Hermitian adjoint.
    % Prefer A.adjoint() so InvLinearOp inputs preserve invertibility.

    properties (SetAccess = private)
        % Constructor validates the type and scalar size.
        % Avoid implicit default construction of an abstract class.
        operator
    end

    methods
        function obj = AjointOfLinearOp(operator)
            arguments
                operator (1,1) LinearOp
            end
            obj@LinearOp(operator.dimRange,operator.dimDomain);
            obj.operator = operator;
        end

        function y = fwd(obj,x)
            obj.validateDomainVector(x);
            y = obj.operator.adj(x);
        end

        function x = adj(obj,y)
            obj.validateRangeVector(y);
            x = obj.operator.fwd(y);
        end

        function A = matrixRep(obj)
            A = obj.operator.matrixRep();
            A = A';
        end

        function adjointOp = adjoint(obj)
            % (A*)* = A; avoid accumulating wrappers.
            adjointOp = obj.operator;
        end
    end
end
