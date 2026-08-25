function [K,B,C,cMax] = sampleMedium1D(xGrid,medium,options)
% sampleMedium1D samples a continuous acoustic medium onto 1D grids.
% Kappa is sampled on the supplied grid, while beta may be sampled on
% either the supplied primal grid or its associated dual grid. Optional
% outputs C and cMax contain the primal-grid wave speed and its maximum.

arguments
    xGrid (1,1) GridSpace1D
    medium (1,1) AcousticMedium1D
    options.BetaGrid (1,1) string ...
        {mustBeMember(options.BetaGrid,["primal","dual"])} = "dual"
end

K = GridFunction(xGrid, ...
    Label=medium.kappaLabel, ...
    Units=medium.kappaUnits);
K.setValues(@(x) medium.evaluateKappa(x));

switch options.BetaGrid
    case "primal"
        betaGrid = xGrid;

    case "dual"
        if ~isa(xGrid.x,"Grid1DUnifPrimal")
            error("sampleMedium1D:PrimalGridRequired", ...
                "A primal uniform grid is required to construct its dual.");
        end

        betaGrid = GridSpace1D(Grid1DUnifDual(xGrid.x));
end

B = GridFunction(betaGrid, ...
    Label=medium.betaLabel, ...
    Units=medium.betaUnits);
B.setValues(@(x) medium.evaluateBeta(x));

if nargout >= 3
    C = GridFunction(xGrid, ...
        Label=medium.waveSpeedLabel, ...
        Units=medium.waveSpeedUnits);
    C.setValues(@(x) medium.evaluateWaveSpeed(x));
end

if nargout >= 4
    cMax = max(C.values,[],'all');
end

end
