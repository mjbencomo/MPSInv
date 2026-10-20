function h = plotAcqGeo(gridFun,sourcePoints,receiverPoints,options)
%PLOTACQGEO Plot a Cartesian field and source/receiver locations.
%   h = plotAcqGeo(gridFun,sourceLocations,receiverLocations)
%   h = plotAcqGeo(...,Parent=ax)
%
% In 1D, sourcePoints and receiverPoints are numeric vectors. In 2D,
% they are numeric M-by-2 matrices with paired [x,y] rows.
% The background GridFunction must be defined on a Cartesian grid.
%
% h.background is the background line or image; h.src_rcv contains the
% source and receiver line handles, in that order.

arguments
    gridFun (1,1) GridFunction
    sourcePoints
    receiverPoints
    options.Parent = []
end

if ~isa(gridFun.grid,'GridSpace1D') && ...
        ~isa(gridFun.grid,'GridSpace2D')
    error('plotAcqGeo:InvalidBackground', ...
        'The background GridFunction must have a Cartesian grid space.');
end

dim = gridFun.dim;
sourceLocations = validateLocations( ...
    sourcePoints,dim,'plotAcqGeo:InvalidSource');
receiverLocations = validateLocations( ...
    receiverPoints,dim,'plotAcqGeo:InvalidReceivers');

if isempty(options.Parent)
    ax = axes(figure);
elseif isgraphics(options.Parent,'axes')
    ax = options.Parent;
else
    error('plotAcqGeo:InvalidParent', ...
        'Parent must be a valid axes object.');
end

h = struct();
hold(ax,'on');
switch dim
    case 1
        yyaxis(ax,'left');
        h.background = plot(ax,gridFun.grid.x.pts,gridFun.values(:),'-k');
        ylabel(ax,gridFun.formattedLabel());
        ax.YAxis(1).Color = h.background.Color;

        yyaxis(ax,'right');
        h.src_rcv = [ ...
            line(ax,'XData',sourceLocations, ...
                'YData',zeros(size(sourceLocations)), ...
                'LineStyle','none','Marker','x','Color','r'), ...
            line(ax,'XData',receiverLocations, ...
                'YData',zeros(size(receiverLocations)), ...
                'LineStyle','none','Marker','o','Color','b')];
        yticks(ax,[]);
        ax.YAxis(2).Color = "k";
        xlabel(ax,'x');

    case 2
        h.background = imagesc(ax,gridFun.grid.x.pts, ...
            gridFun.grid.y.pts,gridFun.values.');
        axis(ax,'xy');
        xlabel(ax,'x');
        ylabel(ax,'y');
        cb = colorbar(ax);
        cb.Label.String = gridFun.formattedLabel();

        h.src_rcv = [ ...
            line(ax,'XData',sourceLocations(:,1), ...
                'YData',sourceLocations(:,2), ...
                'LineStyle','none','Marker','x','Color','r'), ...
            line(ax,'XData',receiverLocations(:,1), ...
                'YData',receiverLocations(:,2), ...
                'LineStyle','none','Marker','o','Color','b')];
        set(h.src_rcv,'MarkerSize',6);
end
hold(ax,'off');
end

function locations = validateLocations(input,dim,errorId)
    if ~isnumeric(input) || ~isreal(input) || any(~isfinite(input),'all')
        error(errorId,'Locations must contain finite real numeric values.');
    end

    if isempty(input)
        locations = zeros(0,dim);
    elseif dim == 1 && isvector(input)
        locations = input(:);
    elseif dim == 2 && ismatrix(input) && size(input,2) == 2
        locations = input;
    else
        error(errorId, ...
            '1D locations must be a vector; 2D locations must be M-by-2.');
    end
end
