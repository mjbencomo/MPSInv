classdef (Abstract) MultipoleComponent
<<<<<<< HEAD
    % MultipoleComponent describes the time-independent part of a
    % multipole source term.
    %
    % A component specifies source geometry and spatial discretization.
    % Concrete 1D and 2D subclasses add a time function to construct the
    % corresponding MultipoleTerm object used by the acoustic solvers.

    properties (SetAccess = private)
        approximationOrder (1,1) double = 4
    end

    methods
        function obj = MultipoleComponent(approximationOrder)
            arguments
                approximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            obj.approximationOrder = approximationOrder;
        end
    end

    methods (Abstract)
        term = withTimeFunction(obj,timeFunction,options)
        % Combine this spatial component with a temporal function.
=======
    % MultipoleComponent describes the time-independent spatial part of a
    % multipole source term.
    %
    % A component contains only continuous source information. Numerical
    % choices, such as the approximation order and computational grid, are
    % supplied when createStencil is called.

    methods (Abstract)
        text = description(obj)
        % Return a readable description of the spatial component.

        [indices,weights,varargout] = createStencil(obj,varargin)
        % Discretize the component on a supplied grid.
>>>>>>> b0d5b1dcea97a019fdcf42ee218d06a571246d9f
    end
end
