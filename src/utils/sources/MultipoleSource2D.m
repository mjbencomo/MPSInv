classdef MultipoleSource2D
    % MultipoleSource2D discretizes sums of 2D multipole terms.
    %
    % A term with derivativeOrder=[sx,sy] is discretized as the tensor
    % product of one-dimensional MPS approximations in x and y. The result
    % is an AcousticSource2D containing localized terms on the appropriate
    % pressure, x-velocity, and y-velocity staggered grids.

    properties (SetAccess = private)
        terms
        approximationOrder (1,1) double
    end

    methods
        function obj = MultipoleSource2D(terms,options)
            arguments
                terms
                options.ApproximationOrder (1,1) double ...
                    {mustBeInteger,mustBePositive} = 4
            end

            if isempty(terms)
                terms = MultipoleTerm2D.empty(1,0);
            elseif iscell(terms)
                if any(~cellfun(@(term) isa(term,'MultipoleTerm2D'),terms))
                    error('MultipoleSource2D:InvalidTerms', ...
                        'Every term must be a MultipoleTerm2D object.');
                end
                terms = [terms{:}];
            elseif ~isa(terms,'MultipoleTerm2D')
                error('MultipoleSource2D:InvalidTerms', ...
                    'terms must contain MultipoleTerm2D objects.');
            end

            obj.terms = reshape(terms,1,[]);
            obj.approximationOrder = options.ApproximationOrder;
        end

        function source = discretize(obj,pressureGrid,velocityXGrid,velocityYGrid)
            % Discretize all terms on their respective staggered grids.
            arguments
                obj
                pressureGrid (1,1) GridSpace2D
                velocityXGrid (1,1) GridSpace2D
                velocityYGrid (1,1) GridSpace2D
            end

            source = AcousticSource2D.zero();
            q = obj.approximationOrder;

            for k = 1:numel(obj.terms)
                term = obj.terms(k);
                switch term.targetField
                    case "pressure"
                        [indices,weights,xIndices,yIndices] = ...
                            obj.tensorProductStencil(pressureGrid,term,q);

                        % The baseline non-PML solver overwrites pressure
                        % along every outer edge. Require fully interior
                        % pressure support, as in MultipoleSource1D.
                        if any(xIndices == 1 | ...
                                xIndices == pressureGrid.N(1)) || ...
                                any(yIndices == 1 | ...
                                yIndices == pressureGrid.N(2))
                            error( ...
                                'MultipoleSource2D:PressureStencilAtBoundary', ...
                                ['The pressure multipole stencil must ' ...
                                 'exclude all pressure-grid edges.']);
                        end

                        source.addLocalizedPressureTerm( ...
                            pressureGrid,indices,weights, ...
                            @(t) term.evaluateTime(t));

                    case "velocityX"
                        [indices,weights] = ...
                            obj.tensorProductStencil(velocityXGrid,term,q);
                        source.addLocalizedVelocityXTerm( ...
                            velocityXGrid,indices,weights, ...
                            @(t) term.evaluateTime(t));

                    case "velocityY"
                        [indices,weights] = ...
                            obj.tensorProductStencil(velocityYGrid,term,q);
                        source.addLocalizedVelocityYTerm( ...
                            velocityYGrid,indices,weights, ...
                            @(t) term.evaluateTime(t));
                end
            end
        end

        function value = numberOfTerms(obj)
            value = numel(obj.terms);
        end
    end

    methods (Access = private)
        function [indices,weights,xIndices,yIndices] = ...
                tensorProductStencil(~,grid,term,q)
            [~,xIndices,xWeights] = MPSappx( ...
                grid.x.pts,term.location(1),q,term.derivativeOrder(1));
            [~,yIndices,yWeights] = MPSappx( ...
                grid.y.pts,term.location(2),q,term.derivativeOrder(2));

            [IX,IY] = ndgrid(xIndices,yIndices);
            weights = xWeights*yWeights.';
            indices = sub2ind(grid.N,IX,IY);
            indices = indices(:);
            weights = weights(:);
        end
    end
end
