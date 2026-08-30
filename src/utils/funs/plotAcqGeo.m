function h= plotAcqGeo(gridFun,source_mps,receivers_gsp,options)
%PLOTACQGEO Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    gridFun (1,1) GridFunction
    source_mps (1,1) 
    receivers_gsp (1,1) GridSpace
    options.Parent = []
end

if gridFun.dim==1
    % 1D case ------------------------------------------------------------%
    if ~isa(source_mps,"MultipoleTerm1D") && ~isa(source_mps,"MultipoleSeries1D")
        error('plotAcqGeo:InvalidSource', ...
            'Source must be of type MultipoleTerm1D or MultipoleSeries1D for 1D case.');
    end
    if ~isa(receivers_gsp, "GridSpace1D")
        error('plotAcqGeo:InvalidReceivers', ...
            'Receiver must be of type GridSpace1D for 1D case.');
    end

    % Extracting source locations
    source_locX = [];
    if isa(source_mps,"MultipoleTerm")
        source_locX = source_mps.location;
    else
        source_locX = source_mps.locations;
    end

    % Extracting receiver locations
    receiver_locX = receivers_gsp.mesh;
    
    % Create or validate figure axis
    if isempty(options.Parent)
        ax = axes(figure);
    elseif isgraphics(options.Parent,'axes')
        ax = options.Parent;
    else
        error('plotAcqGeo:InvalidParent', ...
            'Parent must be a valid axes object.');
    end

    % Plotting background grid function
    yyaxis(ax,"left");
    h.background = plot(ax, ...
        gridFun.grid.mesh, ...
        gridFun.values(:),"-k");
    ylabel(ax,gridFun.formattedLabel);
    ax.YAxis(1).Color = h.background.Color; 

    % Plotting source location
    yyaxis(ax,"right")
    h.src_rcv = plot(ax, ...
        source_locX,0,"xr", ...
        receiver_locX,0,"ob");
    yticks(ax,[]);
    ax.YAxis(2).Color = "k"; 

    xlabel(ax,'x');


elseif gridFun.dim==2
    % 2D case ------------------------------------------------------------%
    if ~isa(souce_mps,"MultipoleTerm2D") && ~isa(source_mps,"MultipoleSeries2D")
        error('plotAcqGeo:InvalidSource', ...
            'Source must be of type MultipoleTerm2D or MultipoleSeries2D for 2D case.');
    end

    if ~isa(receivers_gsp, "GridSpace2D")
        error('plotAcqGeo:InvalidReceivers', ...
            'Receiver must be of type GridSpace2D for 2D case.');
    end

else
    error('plotAcqGeo:InvalidDimension', ...
        'Plotting is only supported for 1D and 2D grids.');
end

end