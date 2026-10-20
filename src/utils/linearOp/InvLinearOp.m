classdef (Abstract) InvLinearOp < LinearOp
    % InvLinearOp is an invertible finite-dimensional linear operator.
    %
    % Since an invertible operator maps a vector space onto itself, its
    % domain and range must have the same dimension.
    %
    % Subclasses must implement:
    %   fwd(obj,x)       - apply the operator
    %   adj(obj,y)       - apply its adjoint
    %   matrixRep(obj)   - return its matrix representation
    %   inverse(obj)     - return an operator representing its inverse

    methods
        function obj = InvLinearOp(dim)
            arguments
                dim (1,1) double ...
                    {mustBeInteger,mustBePositive}
            end

            obj@LinearOp(dim,dim);
        end
    end

    methods (Abstract)
        inverseOp = inverse(obj)
        % Return a LinearOp representing the inverse of this operator.
        %
        % The returned operator should satisfy, up to numerical error,
        %
        %   inverseOp.fwd(obj.fwd(x)) = x
        %   obj.fwd(inverseOp.fwd(x)) = x
    end
end