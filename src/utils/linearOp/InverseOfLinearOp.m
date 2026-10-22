classdef InverseOfLinearOp < InvLinearOp
    % Wrapper representing A^(-1), using A's inverse actions.
    % Construct through A.inverse(). Only matrixRep assembles a matrix.

    properties (SetAccess = private)
        % Constructor validates the type and scalar size.
        % Avoid implicit default construction of an abstract class.
        operator
    end

    methods
        function obj = InverseOfLinearOp(operator)
            arguments
                operator (1,1) InvLinearOp
            end
            obj@InvLinearOp(operator.dimDomain);
            obj.operator = operator;
        end

        function y = fwd(obj,x)
            obj.validateDomainVector(x);
            y = obj.operator.inv(x);
        end

        function x = adj(obj,y)
            obj.validateRangeVector(y);
            x = obj.operator.invAdj(y);
        end

        function x = inv(obj,b)
            obj.validateRangeVector(b);
            x = obj.operator.fwd(b);
        end

        function x = invAdj(obj,b)
            obj.validateDomainVector(b);
            x = obj.operator.adj(b);
        end

        function A = matrixRep(obj)
            % Explicit inverse requested: solve against the identity.
            % No call to inv(A), which would obscure action/operator roles.
            A = obj.operator.matrixRep();
            A = A\eye(obj.dimDomain,'like',A);
        end

        function inverseOp = inverse(obj)
            % (A^(-1))^(-1) = A; do not add another wrapper.
            inverseOp = obj.operator;
        end
    end
end
