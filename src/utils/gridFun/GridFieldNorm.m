classdef GridFieldNorm
    % Norms reduce storage axes, never use a matrix spectral norm.
    % Automatic quadrature: trapezoidal primal/time and midpoint dual.
    % Receiver measures and unresolved singleton widths require user weights.
    methods (Static)
        function result = compute(obj,p,options)
            if ~ismember(p,[1,2,Inf])
                error('GridFieldNorm:InvalidOrder','p must be 1, 2 or Inf.');
            end
            series = isa(obj,'GridFunctionTimeSeries');
            receivers = isa(obj.grid,'ReceiverSpace1D') || ...
                isa(obj.grid,'ReceiverSpace2D');
            if receivers || obj.dim == 1
                spatialSize = obj.numPts;
            else
                spatialSize = obj.N;
            end
            ns = numel(spatialSize);
            shape = spatialSize;
            if series
                shape = [shape,obj.Nt];
            end
            nAxes = numel(shape);
            dims = options.Dimension;
            if ~isempty(dims)
                if options.Over ~= ""
                    error('GridFieldNorm:ConflictingReduction', ...
                        'Specify Over or Dimension, not both.');
                end
                if ~isreal(dims) || ~isvector(dims) || any(~isfinite(dims)) || ...
                        any(dims ~= fix(dims)) || any(dims < 1 | dims > nAxes) || ...
                        numel(unique(dims)) ~= numel(dims)
                    error('GridFieldNorm:InvalidDimension', ...
                        'Dimension must contain distinct valid storage axes.');
                end
                dims = sort(dims(:).');
                scope = "axes " + strjoin(string(dims),",");
            else
                over = options.Over;
                if over == ""
                    over = "all";
                end
                switch over
                    case "all"
                        dims = 1:nAxes;
                    case "space"
                        dims = 1:ns;
                    case "time"
                        if ~series
                            error('GridFieldNorm:NoTimeAxis', ...
                                'A GridFunction has no time axis.');
                        end
                        dims = ns+1;
                    otherwise
                        error('GridFieldNorm:InvalidReduction', ...
                            'Over must be space, time or all.');
                end
                scope = over;
            end
            weighted = options.Weighted || ~isempty(options.Weights);
            nr = prod(shape(dims));
            weights = ones(nr,1);
            if weighted
                weights = GridFieldNorm.weights(obj,options.Weights, ...
                    dims,shape,ns,receivers);
            elseif options.MeasureUnits ~= ""
                error('GridFieldNorm:UnweightedMeasureUnits', ...
                    'MeasureUnits applies only to weighted norms.');
            end

            % Collapse selected axes into rows and all retained axes into
            % columns. Singleton axes retain their original positions.
            retained = setdiff(1:nAxes,dims,'stable');
            dataShape = shape;
            if numel(dataShape) == 1
                dataShape = [dataShape,1];
            end
            permutation = [dims,retained];
            if nAxes == 1
                permutation = [permutation,2];
            end
            values = reshape(double(obj.values),dataShape);
            values = reshape(permute(values,permutation),nr,[]);
            positive = weights > 0;
            values = abs(values(positive,:));
            weights = weights(positive);
            if isinf(p)
                reduced = max(values,[],1,'includenan');
            else
                % Apply pth-root weights before scaling so tiny samples
                % with large weights are not lost during exponentiation.
                values = values .* weights.^(1/p);
                scale = max(values,[],1,'includenan');
                scaled = values ./ scale;
                scaled(:,scale == 0) = 0;
                reduced = scale .* sum(scaled.^p,1).^(1/p);
                reduced(isinf(scale)) = Inf;
            end
            outShape = shape;
            outShape(dims) = 1;
            if numel(outShape) == 1
                outShape = [outShape,1];
            end
            reduced = reshape(reduced,outShape);
            timeReduced = series && ismember(ns+1,dims);
            spaceReduced = all(ismember(1:ns,dims));
            if spaceReduced && (~series || timeReduced)
                result = reduced(1);
                return;
            end

            grid = GridFieldNorm.reducedGrid(obj.grid, ...
                intersect(dims,1:ns),receivers);
            if series && ~timeReduced
                result = GridFunctionTimeSeries(grid,obj.timePoints);
            else
                result = GridFunction(grid);
            end
            result.copyMetadataFrom(obj);
            kind = "Discrete ";
            if weighted
                kind = "Weighted ";
            end
            result.label = kind + scope + " " + string(p) + ...
                "-norm of " + obj.label;
            % Coordinate units are not stored by GridSpace. The caller may
            % supply the combined units of the integration measure.
            if weighted && ~isinf(p) && options.MeasureUnits ~= ""
                factor = "(" + options.MeasureUnits + ")^(1/" + string(p) + ")";
                if obj.units == ""
                    result.units = factor;
                else
                    result.units = "(" + obj.units + ")*" + factor;
                end
            end
            result.setValues(reduced);
        end
    end

    methods (Static, Access = private)
        function w = weights(obj,provided,dims,shape,ns,receivers)
            if isempty(provided)
                axesWeights = cell(1,numel(dims));
                for j = 1:numel(dims)
                    d = dims(j);
                    if d > ns
                        axesWeights{j} = GridFieldNorm.trapezoidal(obj.timePoints);
                    elseif receivers
                        error('GridFieldNorm:ReceiverWeightsRequired', ...
                            'Spatial receiver norms require explicit quadrature weights.');
                    else
                        coordinate = obj.grid.x;
                        if d == 2
                            coordinate = obj.grid.y;
                        end
                        if isa(coordinate,'Grid1DUnifDual')
                            axesWeights{j} = coordinate.h*ones(coordinate.N,1);
                        else
                            axesWeights{j} = GridFieldNorm.trapezoidal(coordinate.pts);
                        end
                    end
                end
                provided = axesWeights;
            end
            if iscell(provided)
                if numel(provided) ~= numel(dims)
                    error('GridFieldNorm:InvalidWeights', ...
                        'Supply one weight vector per reduced axis in ascending order.');
                end
                w = 1;
                for j = 1:numel(dims)
                    v = provided{j};
                    GridFieldNorm.validateWeights(v);
                    if ~isvector(v) || numel(v) ~= shape(dims(j))
                        error('GridFieldNorm:InvalidWeights', ...
                            'Each weight vector must match its reduced axis.');
                    end
                    w = kron(double(v(:)),w);
                end
            else
                GridFieldNorm.validateWeights(provided);
                expected = shape(dims);
                if numel(expected) == 1
                    valid = isvector(provided) && numel(provided) == expected;
                else
                    % MATLAB suppresses trailing singleton dimensions.
                    expected(end+1) = 1;
                    while numel(expected) > 2 && expected(end) == 1
                        expected(end) = [];
                    end
                    valid = isequal(size(provided),expected);
                end
                if ~valid
                    error('GridFieldNorm:InvalidWeights', ...
                        'Weight array size must match the reduced axes in ascending order.');
                end
                w = double(provided(:));
            end
            GridFieldNorm.validateWeights(w);
        end

        function validateWeights(w)
            if ~isnumeric(w) || isempty(w) || ~isreal(w) || ...
                    any(~isfinite(w(:))) || any(w(:) < 0) || ~any(w(:) > 0)
                error('GridFieldNorm:InvalidWeights', ...
                    'Weights must be finite, real, nonnegative and have positive total weight.');
            end
        end

        function w = trapezoidal(points)
            points = points(:);
            if numel(points) < 2
                error('GridFieldNorm:SingletonWeightsRequired', ...
                    'An integration axis with one sample requires an explicit weight.');
            end
            h = diff(points);
            w = [h(1);h(1:end-1)+h(2:end);h(end)]/2;
        end

        function grid = reducedGrid(grid,dims,receivers)
            if isempty(dims)
                return;
            end
            if receivers
                % Representative centroid only, not an observation location.
                if isa(grid,'ReceiverSpace1D')
                    grid = ReceiverSpace1D(mean(grid.points));
                else
                    grid = ReceiverSpace2D(mean(grid.points,1));
                end
                return;
            end
            x = grid.x;
            if ismember(1,dims)
                x = Grid1D(x.dom,mean(x.dom));
            end
            if grid.dim == 1
                grid = GridSpace1D(x);
            else
                y = grid.y;
                if ismember(2,dims)
                    y = Grid1D(y.dom,mean(y.dom));
                end
                grid = GridSpace2D(x,y);
            end
        end
    end
end
