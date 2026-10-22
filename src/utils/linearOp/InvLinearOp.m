classdef (Abstract) InvLinearOp < LinearOp
    % Invertible finite-dimensional operator with inverse actions.
    %
    % Subclasses implement fwd, adj, matrixRep, inv, and invAdj.
    %   inv(b)    solves A*x = b.
    %   invAdj(b) solves A'*x = b (Hermitian adjoint).
    %
    % inverse() and adjoint() return InvLinearOp wrappers without
    % constructing matrices or performing solves at wrapper creation.
    % Domain and range have equal dimension but may represent different
    % coordinates. These are inverse actions, not pseudoinverse actions.

    methods
        function obj = InvLinearOp(dim)
            arguments
                dim (1,1) double ...
                    {mustBeInteger,mustBePositive,mustBeFinite}
            end
            obj@LinearOp(dim,dim);
        end

        function inverseOp = inverse(obj)
            inverseOp = InverseOfLinearOp(obj);
        end

        function adjointOp = adjoint(obj)
            adjointOp = AdointInvOfLinearOp(obj);
        end
    end

    methods (Abstract)
        x = inv(obj,b)
        % Apply A^(-1). b lies in A's range; x lies in A's domain.

        x = invAdj(obj,b)
        % Apply (A')^(-1). b lies in A's domain; x lies in A's range.
    end
end
