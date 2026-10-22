classdef AdointInvOfLinearOp < InvLinearOp
    % Wrapper representing A', preserving invertibility.
    % Construction and actions do not assemble matrices.

    properties (SetAccess = private)
        % Constructor validates the type and scalar size.
        % Avoid implicit default construction of an abstract class.
        operator
    end

    methods
        function obj = AdointInvOfLinearOp(operator)
            arguments
                operator (1,1) InvLinearOp
            end
            obj@InvLinearOp(operator.dimDomain);
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

        function x = inv(obj,b)
            obj.validateRangeVector(b);
            x = obj.operator.invAdj(b);
        end

        function x = invAdj(obj,b)
            obj.validateDomainVector(b);
            x = obj.operator.inv(b);
        end

        function A = matrixRep(obj)
            A = obj.operator.matrixRep();
            A = A';
        end

        function adjointOp = adjoint(obj)
            % (A')' = A. inverse() is inherited from InvLinearOp.
            adjointOp = obj.operator;
        end
    end
end
