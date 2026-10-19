classdef MatrixFreeInvOpForTest < InvLinearOp
    % Test double: diagonal actions are available; matrix assembly fails.
    properties (SetAccess = private)
        diagonal (:,1) double
    end
    methods
        function obj = MatrixFreeInvOpForTest(diagonal)
            obj@InvLinearOp(numel(diagonal));
            obj.diagonal = diagonal(:);
        end
        function y = fwd(obj,x)
            obj.validateDomainVector(x);
            y = obj.diagonal.*x;
        end
        function x = adj(obj,y)
            obj.validateRangeVector(y);
            x = conj(obj.diagonal).*y;
        end
        function x = inv(obj,b)
            obj.validateRangeVector(b);
            x = b./obj.diagonal;
        end
        function x = invAdj(obj,b)
            obj.validateDomainVector(b);
            x = b./conj(obj.diagonal);
        end
        function A = matrixRep(~) %#ok<STOUT>
            error('MatrixFreeInvOpForTest:MatrixForbidden', ...
                'This test operator intentionally forbids matrix assembly.');
        end
    end
end
