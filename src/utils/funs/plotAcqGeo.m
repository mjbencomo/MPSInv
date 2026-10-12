function h = plotAcqGeo(gridFun, source_mps, receivers_gsp, options)
%PLOTACQGEO Plot a Cartesian-grid background and acquisition geometry.
%   h = plotAcqGeo(gridFun, source_mps, receivers_gsp)
%   h = plotAcqGeo(..., Parent=ax)
%
%   h is a struct with fields:
%     - background : handle to background plot (line or image)
%     - src_rcv    : two line handles: sources and receivers

arguments
    gridFun (1,1) GridFunction
    source_mps
    receivers_gsp (1,1) GridSpace
    options.Parent = []
end

h = struct();

if ~isa(gridFun.grid,'GridSpace1D') && ...
        ~isa(gridFun.grid,'GridSpace2D')
    error('plotAcqGeo:InvalidBackground', ...
        'The background GridFunction must have a Cartesian grid space.');
end

% Create or validate axis
if isempty(options.Parent)
    ax = axes(figure);
elseif isgraphics(options.Parent,'axes')
    ax = options.Parent;
else
    error('plotAcqGeo:InvalidParent', 'Parent must be a valid axes object.');
end
hold(ax,'on')

switch gridFun.dim
    case 1
        % Validate types
        if ~isa(source_mps,'MultipoleSrcTerm1D') && ~isa(source_mps,'MultipoleSrcSeries1D')
            error('plotAcqGeo:InvalidSource', ...
                'Source must be MultipoleTerm1D or MultipoleSeries1D for 1D.');
        end
        if ~isa(receivers_gsp,'GridSpace1D') && ...
                ~isa(receivers_gsp,'ReceiverSpace1D')
            error('plotAcqGeo:InvalidReceivers', ...
                ['Receivers must be a GridSpace1D or ' ...
                 'ReceiverSpace1D for 1D plotting.']);
        end

        % Source locations (column vector)
        if isa(source_mps,'MultipoleSrcTerm1D')
            src_locX = source_mps.location;
        else
            src_locX = source_mps.locations;
        end
        src_locX = src_locX(:);

        if isa(receivers_gsp,'ReceiverSpace1D')
            rcv_locX = receivers_gsp.points;
        else
            rcv_locX = receivers_gsp.mesh();
            rcv_locX = rcv_locX(:);
        end

        % Plot background (left axis)
        yyaxis(ax,'left');
        h.background = plot(ax, gridFun.grid.x.pts, gridFun.values(:), '-k');
        ylabel(ax, gridFun.formattedLabel());
        ax.YAxis(1).Color = h.background.Color;

        % Plot sources and receivers on right axis (y=0)
        yyaxis(ax,'right');
        h.src_rcv = plot(ax, src_locX, zeros(size(src_locX)), 'xr', ...
                             rcv_locX, zeros(size(rcv_locX)), 'ob');
        yticks(ax,[]);
        ax.YAxis(2).Color = "k";
        xlabel(ax,'x');

    case 2
        % Validate types
        if ~isa(source_mps,'MultipoleSrcTerm2D') && ~isa(source_mps,'MultipoleSrcSeries2D')
            error('plotAcqGeo:InvalidSource', ...
                'Source must be MultipoleTerm2D or MultipoleSeries2D for 2D.');
        end
        if ~isa(receivers_gsp,'GridSpace2D') && ...
                ~isa(receivers_gsp,'ReceiverSpace2D')
            error('plotAcqGeo:InvalidReceivers', ...
                ['Receivers must be a GridSpace2D or ' ...
                 'ReceiverSpace2D for 2D plotting.']);
        end

        % Source locations (Nx2)
        if isa(source_mps,'MultipoleSrcTerm2D')
            source_loc = source_mps.location;
        else
            source_loc = source_mps.locations;
        end
        if isempty(source_loc)
            src_locX = zeros(0,1); src_locY = zeros(0,1);
        else
            src_locX = source_loc(:,1);
            src_locY = source_loc(:,2);
        end

        if isa(receivers_gsp,'ReceiverSpace2D')
            % Keep each paired [x,y] row as one receiver marker.
            rcv_locX = receivers_gsp.points(:,1);
            rcv_locY = receivers_gsp.points(:,2);
        else
            [rcv_locX,rcv_locY] = receivers_gsp.mesh();
            rcv_locX = rcv_locX(:);
            rcv_locY = rcv_locY(:);
        end

        % Plot background image
        h.background = imagesc(ax, gridFun.grid.x.pts, gridFun.grid.y.pts, gridFun.values.');
        axis(ax,'xy');
        xlabel(ax,'x');
        ylabel(ax,'y');
        cb = colorbar(ax);
        cb.Label.String = gridFun.formattedLabel();

        % Overlay sources and receivers
        h.src_rcv = plot(ax, src_locX, src_locY, 'xr', rcv_locX, rcv_locY, 'ob');
        set(h.src_rcv, 'MarkerSize', 6);

    otherwise
        error('plotAcqGeo:InvalidDimension', 'Only 1D and 2D grids are supported.');
end

hold(ax,'off')
end
