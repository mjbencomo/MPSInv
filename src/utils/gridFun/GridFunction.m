% GridFunction represents numerical values defined on a spatial grid.
% Values are initialized to zero and may be set using data or a function.

classdef GridFunction < handle
    properties (SetAccess = protected)
        grid
        values
    end

    properties
        label (1,1) string = "Value"
        units (1,1) string = ""
    end

    properties (Dependent)
        N
        numPts
        dim
    end

    %%
    methods
        
        %%%%%%%%
        function obj = GridFunction(grid,options)
            % Construct a zero-valued, optionally labeled grid function.

            arguments
                grid (1,1) GridSpace
                options.Label (1,1) string = "Value"
                options.Units (1,1) string = ""
            end

            obj.grid = grid;
            obj.label = options.Label;
            obj.units = options.Units;

            % Initialize values
            if isa(grid,'ReceiverSpace1D') || isa(grid,'ReceiverSpace2D')
                obj.values = zeros(grid.numPts,1);

            elseif grid.dim == 1
                obj.values = zeros(grid.N,1);

            elseif grid.dim == 2
                obj.values = zeros(grid.N);

            else
                error('GridFunction:InvalidDimension', ...
                    'Unsupported spatial dimension.');
            end
        end

        %%%%%%%%
        function setValues(obj,data)
            % Set grid-function values using an array or function handle.

            arguments
                obj
                data
            end

            if isa(data,'function_handle')
                obj.setValuesFromFunction(data);

            elseif isnumeric(data)
                obj.setValuesFromArray(data);

            else
                error('GridFunction:InvalidData', ...
                    'Input must be a numeric array or function handle.');
            end
        end

        %%%%%%%%
        function value = get.N(obj)
            if isa(obj.grid,'ReceiverSpace1D') || ...
                    isa(obj.grid,'ReceiverSpace2D')
                error('GridFunction:NUndefinedForReceivers', ...
                    'Receiver spaces have numPts, but no Cartesian grid size N.');
            end
            value = obj.grid.N;
        end

        function value = get.numPts(obj)
            value = obj.grid.numPts;
        end

        %%%%%%%%
        function value = get.dim(obj)
            value = obj.grid.dim;
        end

        %%%%%%%%
        function value = formattedLabel(obj)
            % Return the field label with units when units are provided.

            if strlength(obj.units) == 0
                value = obj.label;
            else
                value = obj.label + " (" + obj.units + ")";
            end
        end

        %%%%%%%%
        function h = plot(obj,varargin)
            if isa(obj.grid,'ReceiverSpace1D')
                h = plot(obj.grid.points,obj.values,'o',varargin{:});
                xlabel('x');
                ylabel(obj.formattedLabel());

            elseif isa(obj.grid,'ReceiverSpace2D')
                h = scatter(obj.grid.points(:,1),obj.grid.points(:,2), ...
                    [],obj.values,'filled',varargin{:});
                xlabel('x');
                ylabel('y');
                cb = colorbar;
                cb.Label.String = obj.formattedLabel();

            elseif obj.dim == 1
                if obj.grid.N == 1
                    h = plot(obj.grid.x.pts,obj.values, ...
                        'o',varargin{:});
                else
                    h = plot(obj.grid.x.pts,obj.values, ...
                        varargin{:});
                end
                xlabel('x');
                ylabel(obj.formattedLabel());

            elseif obj.dim == 2
                % [X,Y] = obj.grid.mesh();
                h = imagesc( ...
                        obj.grid.x.pts, ...
                        obj.grid.y.pts, ...
                        obj.values.', ...
                        varargin{:});
                xlabel('x');
                ylabel('y');
                cb = colorbar;
                cb.Label.String = obj.formattedLabel();

            else
                error('GridFunction:InvalidDimension', ...
                    'Plotting is only supported for 1D and 2D grids.');
            end
        end

        %%%%%%%%
        function interpolateFrom(obj,source,method)
            % Interpolate values from another GridFunction onto this grid.
            %
            %   obj.interpolateFrom(source)
            %   obj.interpolateFrom(source,method)
            %
            % The source and target GridFunctions must have the same spatial
            % dimension. The default interpolation method is linear. Only
            % numerical values are transferred: obj.label and obj.units
            % remain unchanged because the target metadata is authoritative.

            arguments
                obj
                source (1,1) GridFunction
                method (1,1) string = "linear"
            end

            % Check spatial dimensions
            if obj.dim ~= source.dim
                error('GridFunction:DimensionMismatch', ...
                    'Source and target GridFunctions must have the same dimension.');
            end

            if isa(obj.grid,'ReceiverSpace1D') || ...
                    isa(obj.grid,'ReceiverSpace2D') || ...
                    isa(source.grid,'ReceiverSpace1D') || ...
                    isa(source.grid,'ReceiverSpace2D')
                if isequal(obj.grid,source.grid)
                    obj.setValues(source.values);
                    return;
                end
                error('GridFunction:ReceiverInterpolationUnsupported', ...
                    ['Interpolation involving different receiver or ' ...
                     'Cartesian spaces is not implemented yet.']);
            end

            if obj.dim == 1

                % Target coordinates
                X = obj.grid.x.pts;

                % Interpolate
                data = interp1( ...
                    source.grid.x.pts, ...
                    source.values, ...
                    X, ... 
                    method);

            elseif obj.dim == 2

                % Target coordinates
                [X,Y] = obj.grid.mesh();

                % Interpolating data
                data = interp2( ...
                    source.grid.y.pts, ...
                    source.grid.x.pts, ...
                    source.values, ...
                    Y,X, ...
                    "linear");

            else
                error('GridFunction:InvalidDimension', ...
                    'Interpolation is only supported for 1D and 2D grids.');
            end

            obj.setValues(data);
        end

    end

    %% Arithmetic methods
    % Arithmetic returns independent objects. Binary field operands must
    % have the same class, grid and (for series) stored times. Numeric
    % operands must be scalars; use .* and ./ for pointwise operations.
    methods
        function result = plus(a,b)
            result = GridFieldArithmetic.binary(a,b,"plus");
        end

        function result = minus(a,b)
            result = GridFieldArithmetic.binary(a,b,"minus");
        end

        function result = times(a,b)
            result = GridFieldArithmetic.binary(a,b,"times");
        end

        function result = rdivide(a,b)
            result = GridFieldArithmetic.binary(a,b,"rdivide");
        end

        function result = mtimes(a,b)
            result = GridFieldArithmetic.binary(a,b,"mtimes");
        end

        function result = uminus(a)
            result = GridFieldArithmetic.binary(-1,a,"mtimes");
        end
    end

    %% Private methods
    methods (Access = private)

        %%%%%%%%
        function setValuesFromArray(obj,data)

            if isa(obj.grid,'ReceiverSpace1D') || ...
                    isa(obj.grid,'ReceiverSpace2D')
                if ~isvector(data) || numel(data) ~= obj.numPts
                    error('GridFunction:InvalidSize', ...
                        'Receiver data must contain one value per receiver.');
                end
                obj.values = data(:);

            elseif obj.dim == 1
                if ~isvector(data) || numel(data) ~= obj.N
                    error('GridFunction:InvalidSize', ...
                        'Array must contain one value per grid point.');
                end

                obj.values = data(:);

            elseif obj.dim == 2
                if ~isequal(size(data),obj.N)
                    error('GridFunction:InvalidSize', ...
                        'Array size must match the spatial grid.');
                end

                obj.values = data;
            end
        end

        %%%%%%%%
        function setValuesFromFunction(obj,f)

            if isa(obj.grid,'ReceiverSpace1D')
                data = f(obj.grid.points);

            elseif isa(obj.grid,'ReceiverSpace2D')
                xy = obj.grid.points;
                data = f(xy(:,1),xy(:,2));

            elseif obj.dim == 1
                x = obj.grid.x.pts;

                data = f(x);

            elseif obj.dim == 2
                [X,Y] = obj.grid.mesh();

                data = f(X,Y);
            end

            % Use the array routine for validation
            obj.setValuesFromArray(data);
        end
    end
end
