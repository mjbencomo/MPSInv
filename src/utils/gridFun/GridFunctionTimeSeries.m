classdef GridFunctionTimeSeries < handle
    % GridFunctionTimeSeries stores spatial grid functions at specified times.
    % Values are defined on a fixed spatial grid and a numeric array of
    % output times, which may differ from the computational solver grids.

    properties (SetAccess = protected)
        grid
        times
        values
    end

    properties
        label (1,1) string = "Value"
        units (1,1) string = ""
        timeLabel (1,1) string = "Time"
        timeUnits (1,1) string = ""
    end

    properties (Dependent)
        N
        Nt
        dim
    end

    %% Public methods
    methods
        %%%%%%%%
        function obj = GridFunctionTimeSeries(grid,times,options)
            % Construct a zero-valued, optionally labeled time series.

            arguments
                grid (1,1) GridSpace
                times double {mustBeFinite}
                options.Label (1,1) string = "Value"
                options.Units (1,1) string = ""
                options.TimeLabel (1,1) string = "Time"
                options.TimeUnits (1,1) string = ""
            end

            if isempty(times)
                error('GridFunctionTimeSeries:EmptyTimeGrid', ...
                    'At least one output time is required.');
            end

            if ~isvector(times)
                error('GridFunctionTimeSeries:InvalidTimes', ...
                    'Output times must be provided as a vector.');
            end

            times = times(:).';

            if any(diff(times) <= 0)
                error('GridFunctionTimeSeries:InvalidTimes', ...
                    'Output times must be strictly increasing.');
            end

            obj.grid = grid;
            obj.times = times;
            obj.label = options.Label;
            obj.units = options.Units;
            obj.timeLabel = options.TimeLabel;
            obj.timeUnits = options.TimeUnits;

            if grid.dim == 1
                obj.values = zeros(grid.N,numel(times));

            elseif grid.dim == 2
                obj.values = zeros([grid.N,numel(times)]);

            else
                error('GridFunctionTimeSeries:InvalidDimension', ...
                    'Only 1D and 2D spatial grids are supported.');
            end
        end

        %%%%%%%%
        function setSnapshot(obj,timeIndex,data)
            % Store data at one output time.
            % Numeric data must match the spatial grid. A GridFunction is
            % interpolated onto the time-series grid when necessary.

            arguments
                obj
                timeIndex (1,1) double {mustBeInteger,mustBePositive}
                data
            end

            obj.validateTimeIndex(timeIndex);

            if isa(data,'GridFunction')
                if data.dim ~= obj.dim
                    error('GridFunctionTimeSeries:DimensionMismatch', ...
                        ['The source GridFunction and time series must ' ...
                         'have the same spatial dimension.']);
                end

                if isequal(data.grid,obj.grid)
                    data = data.values;
                else
                    target = GridFunction(obj.grid);
                    target.interpolateFrom(data);
                    data = target.values;
                end

            elseif ~isnumeric(data)
                error('GridFunctionTimeSeries:InvalidData', ...
                    'Snapshot data must be a numeric array or GridFunction.');
            end

            if obj.dim == 1
                if ~isvector(data) || numel(data) ~= obj.N
                    error('GridFunctionTimeSeries:InvalidSize', ...
                        'A 1D snapshot must contain one value per grid point.');
                end

                obj.values(:,timeIndex) = data(:);

            elseif obj.dim == 2
                if ~isequal(size(data),obj.N)
                    error('GridFunctionTimeSeries:InvalidSize', ...
                        'A 2D snapshot size must match the spatial grid.');
                end

                obj.values(:,:,timeIndex) = data;
            end
        end

        %%%%%%%%
        function u = snapshot(obj,timeIndex)
            % Return one stored time level as a GridFunction.

            arguments
                obj
                timeIndex (1,1) double {mustBeInteger,mustBePositive}
            end

            obj.validateTimeIndex(timeIndex);
            u = GridFunction(obj.grid, ...
                Label=obj.label, ...
                Units=obj.units);

            if obj.dim == 1
                u.setValues(obj.values(:,timeIndex));
            else
                u.setValues(obj.values(:,:,timeIndex));
            end
        end

        %%%%%%%%
        function h = animate(obj,options)
            % Animate the stored grid function over its output times.
            %
            %   h = U.animate()
            %   h = U.animate(FramePause=0.05,Parent=ax)
            %
            % A 1D series is displayed with plot, while a 2D series is
            % displayed with imagesc. The same graphics object is updated
            % for every frame. Limits are fixed using the complete data set.

            arguments
                obj
                options.FramePause (1,1) double ...
                    {mustBeNonnegative,mustBeFinite} = 0.05
                options.Parent = []
            end

            % Extracting figure axis
            if isempty(options.Parent)
                ax = axes(figure);
            elseif isgraphics(options.Parent,'axes')
                ax = options.Parent;
            else
                error('GridFunctionTimeSeries:InvalidParent', ...
                    'Parent must be a valid axes object.');
            end

            if obj.dim == 1
                if obj.grid.N == 1
                    h = plot(ax,obj.grid.x.pts,obj.values(:,1),'o');
                else
                    h = plot(ax,obj.grid.x.pts,obj.values(:,1));
                end
                xlabel(ax,'x');
                ylabel(ax,obj.formattedLabel());
                ylim(ax,obj.expandedLimits(obj.values));

            elseif obj.dim == 2
                h = imagesc(ax, ...
                    obj.grid.x.pts, ...
                    obj.grid.y.pts, ...
                    obj.values(:,:,1).');
                axis(ax,'xy');
                xlabel(ax,'x');
                ylabel(ax,'y');
                cb = colorbar(ax);
                cb.Label.String = obj.formattedLabel();
                clim(ax,obj.expandedLimits(obj.values));

            else
                error('GridFunctionTimeSeries:InvalidDimension', ...
                    'Animation is only supported for 1D and 2D grids.');
            end

            % Time-loop for animation
            for timeIndex = 1:obj.Nt
                if obj.dim == 1
                    h.YData = obj.values(:,timeIndex);
                else
                    h.CData = obj.values(:,:,timeIndex).';
                end

                title(ax,obj.formattedTime(timeIndex));
                drawnow;

                if options.FramePause > 0 && timeIndex < obj.Nt
                    pause(options.FramePause);
                end
            end
        end

        %%%%%%%%
        function h = plotTimeSeries(obj,options)

            arguments
                obj 
                options.Parent = [] 
            end
        
            % Extracting figure axis
            if isempty(options.Parent)
                ax = axes(figure);
            elseif isgraphics(options.Parent,'axes')
                ax = options.Parent;
            else
                error('GridFunctionTimeSeries:InvalidParent', ...
                    'Parent must be a valid axes object.');
            end

            if obj.dim == 1
                if obj.N == 1
                    h = plot(ax,obj.times,obj.values(1,:));
                    xlabel(ax,obj.timeLabel+" ("+obj.timeUnits+")");
                    ylabel(obj.formattedLabel());
                else
                    h = imagesc(ax, ...
                        1:obj.N, ...
                        obj.times, ...
                        obj.values.');
                    xlabel(ax,"Receiver Index");
                    ylabel(ax,obj.timeLabel+" ("+obj.timeUnits+")");
                    cb = colorbar(ax);
                    cb.Label.String = obj.formattedLabel();
                    clim(ax,obj.expandedLimits(obj.values));
                end

            elseif obj.dim == 2
                % In 2D receivers are numbered first by x, then by y
                % coordinate.
                if obj.N == 1
                    h = plot(ax,obj.times,obj.values(:));
                    xlabel(ax,obj.timeLabel+" ("+obj.timeUnits+")");
                    ylabel(obj.formattedLabel());
                else 
                   Nxy = prod(obj.N);
                   vals = reshape(obj.values,[Nxy,obj.Nt]);
                   h = imagesc(1:Nxy,obj.times,vals.');
                   xlabel(ax,"Receiver Index");
                   ylabel(ax,obj.timeLabel+" ("+obj.timeUnits+")");
                   cb = colorbar(ax);
                   cb.Label.String = obj.formattedLabel();
                   clim(ax,obj.expandedLimits(obj.values));
                end
            end    
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
        function value = formattedTime(obj,timeIndex)
            % Return a formatted label for one stored output time.

            arguments
                obj
                timeIndex (1,1) double {mustBeInteger,mustBePositive}
            end

            obj.validateTimeIndex(timeIndex);
            value = obj.timeLabel + " = " + ...
                sprintf('%.6g',obj.times(timeIndex));

            if strlength(obj.timeUnits) > 0
                value = value + " " + obj.timeUnits;
            end
        end

        %%%%%%%%
        function value = get.N(obj)
            value = obj.grid.N;
        end

        %%%%%%%%
        function value = get.Nt(obj)
            value = numel(obj.times);
        end

        %%%%%%%%
        function value = get.dim(obj)
            value = obj.grid.dim;
        end
    end

    %% Private methods
    methods (Access = private)
        %%%%%%%%
        function limits = expandedLimits(~,data)
            % This function determines range of data values for setting
            % plot limits

            lowerLimit = min(data,[],'all');
            upperLimit = max(data,[],'all');

            if lowerLimit == upperLimit
                padding = 0.05*max(1,abs(lowerLimit));
                limits = [lowerLimit-padding,upperLimit+padding];
            else
                limits = [lowerLimit,upperLimit];
            end
        end

        %%%%%%%%
        function validateTimeIndex(obj,timeIndex)
            % This function checks time index is valid.
            
            if timeIndex > obj.Nt
                error('GridFunctionTimeSeries:InvalidTimeIndex', ...
                    'Time index exceeds the number of stored times.');
            end
        end
    end
end
