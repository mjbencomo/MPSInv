classdef (Abstract) MultipoleComponent
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
    end
end
