function h = plotAcqGeo(gridFun, source_mps, receivers_gsp, options)
%PLOTACQGEO Plot grid function background and acquisition geometry.
%   h = plotAcqGeo(gridFun, source_mps, receivers_gsp)
%   h = plotAcqGeo(..., Parent=ax)
%
%   h is a struct with fields:
%     - background : handle to background plot (line or image)
%     - src_rcv    : handle to source/receiver plot (line object)

arguments
    gridFun (1,1) GridFunction
    source_mps
    receivers_gsp (1,1) GridSpace
    options.Parent = []
end

h = struct();

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
        if ~isa(source_mps,'MultipoleTerm1D') && ~isa(source_mps,'MultipoleSeries1D')
            error('plotAcqGeo:InvalidSource', ...
                'Source must be MultipoleTerm1D or MultipoleSeries1D for 1D.');
        end
        if ~isa(receivers_gsp,'GridSpace1D')
            error('plotAcqGeo:InvalidReceivers', ...
                'Receivers must be a GridSpace1D for 1D plotting.');
        end

        % Source locations (column vector)
        if isa(source_mps,'MultipoleTerm1D')
            src_locX = source_mps.location;
        else
            src_locX = source_mps.locations;
        end
        src_locX = src_locX(:);

        % Receiver locations via mesh()
        rcv_locX = receivers_gsp.mesh();
        rcv_locX = rcv_locX(:);

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
        if ~isa(source_mps,'MultipoleTerm2D') && ~isa(source_mps,'MultipoleSeries2D')
            error('plotAcqGeo:InvalidSource', ...
                'Source must be MultipoleTerm2D or MultipoleSeries2D for 2D.');
        end
        if ~isa(receivers_gsp,'GridSpace2D')
            error('plotAcqGeo:InvalidReceivers', ...
                'Receivers must be a GridSpace2D for 2D plotting.');
        end

        % Source locations (Nx2)
        if isa(source_mps,'MultipoleTerm2D')
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

        % Receiver locations via mesh()
        [rcv_locX, rcv_locY] = receivers_gsp.mesh();
        % GridSpace2D.mesh returns ndgrid(x.pts,y.pts) so rcv_locX,rcv_locY
        % have shape [Nx,Ny]. Flatten for plotting markers.
        rcv_locX = rcv_locX(:);
        rcv_locY = rcv_locY(:);

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