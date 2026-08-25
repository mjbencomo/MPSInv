function [K,Bx,By,C,cMax] = sampleMedium2D(grid,medium,options)
% sampleMedium2D samples a continuous acoustic medium onto 2D grids.
% Kappa is sampled on the supplied primal grid. Beta is sampled on the
% x- and y-staggered grids by default. Optional outputs C and cMax contain
% the primal-grid wave speed and its maximum value.

arguments
    grid (1,1) GridSpace2D
    medium (1,1) AcousticMedium2D
    options.BetaGrid (1,1) string ...
        {mustBeMember(options.BetaGrid,["primal","staggered"])} = ...
        "staggered"
end

K = GridFunction(grid, ...
    Label=medium.kappaLabel, ...
    Units=medium.kappaUnits);
K.setValues(@(x,y) medium.evaluateKappa(x,y));

switch options.BetaGrid
    case "primal"
        bxGrid = grid;
        byGrid = grid;

    case "staggered"
        if ~isa(grid.x,'Grid1DUnifPrimal') || ...
                ~isa(grid.y,'Grid1DUnifPrimal')
            error('sampleMedium2D:PrimalGridRequired', ...
                ['A GridSpace2D containing primal uniform x- and y-grids ' ...
                 'is required to construct the staggered beta grids.']);
        end

        xDual = Grid1DUnifDual(grid.x);
        yDual = Grid1DUnifDual(grid.y);
        bxGrid = GridSpace2D(xDual,grid.y);
        byGrid = GridSpace2D(grid.x,yDual);
end

Bx = GridFunction(bxGrid, ...
    Label=medium.betaLabel, ...
    Units=medium.betaUnits);
By = GridFunction(byGrid, ...
    Label=medium.betaLabel, ...
    Units=medium.betaUnits);
Bx.setValues(@(x,y) medium.evaluateBeta(x,y));
By.setValues(@(x,y) medium.evaluateBeta(x,y));

if nargout >= 4
    C = GridFunction(grid, ...
        Label=medium.waveSpeedLabel, ...
        Units=medium.waveSpeedUnits);
    C.setValues(@(x,y) medium.evaluateWaveSpeed(x,y));
end

if nargout >= 5
    cMax = max(C.values,[],'all');
end

end
