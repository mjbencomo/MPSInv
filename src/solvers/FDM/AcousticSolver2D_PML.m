classdef AcousticSolver2D_PML < AcousticSolver
    % AcousticSolver2D_PML solves 2D acoustics with four PML layers.
    % Pressure is split into x- and y-damped components internally.
    % SpatialOrder may be 2 or 4; leapfrog time stepping is second order.

    properties (SetAccess = protected)
        medium
        pml
        spatialOrder
        extendedGrid
        extendedVelocityXGrid
        extendedVelocityYGrid
        physicalPressureXIndices
        physicalPressureYIndices
        physicalVelocityXXIndices
        physicalVelocityXYIndices
        physicalVelocityYXIndices
        physicalVelocityYYIndices
        K
        Bx
        By
        cMax
        cfl
        sigmaPX
        sigmaPY
        sigmaVX
        sigmaVY
        velocityExtrapolation
    end

    methods
        function obj = AcousticSolver2D_PML( ...
                computationalGrid,computationalTimes, ...
                outputGrid,outputTimes,medium,pml,options)
            arguments
                computationalGrid (1,1) GridSpace2D
                computationalTimes double {mustBeFinite}
                outputGrid (1,1) GridSpace2D
                outputTimes double {mustBeFinite}
                medium (1,1) AcousticMedium2D
                pml (1,1) PML2D
                options.SpatialOrder (1,1) double ...
                    {mustBeMember(options.SpatialOrder,[2,4])} = 2
                options.EnforceCFL (1,1) logical = true
                options.VelocityExtrapolation (1,1) string ...
                    {mustBeMember(options.VelocityExtrapolation, ...
                    ["linear","none"])} = "linear"
            end

            if ~isa(computationalGrid.x,'Grid1DUnifPrimal') || ...
                    ~isa(computationalGrid.y,'Grid1DUnifPrimal')
                error('AcousticSolver2D_PML:PrimalGridRequired', ...
                    ['The physical pressure grid must contain uniform ' ...
                     'primal grids in both directions.']);
            end
            if options.SpatialOrder == 4 && any(computationalGrid.N < 5)
                error('AcousticSolver2D_PML:GridTooSmall', ...
                    ['At least five physical pressure points in each ' ...
                     'direction are required for the fourth-order stencil.']);
            end
            if any(outputGrid.dom(:,1) < computationalGrid.dom(:,1)) || ...
                    any(outputGrid.dom(:,2) > computationalGrid.dom(:,2))
                error('AcousticSolver2D_PML:OutputGridOutsideDomain', ...
                    'Output points must lie in the physical domain.');
            end

            obj@AcousticSolver(computationalGrid,computationalTimes, ...
                outputGrid,outputTimes);

            obj.medium = medium;
            obj.pml = pml;
            obj.spatialOrder = options.SpatialOrder;
            obj.velocityExtrapolation = options.VelocityExtrapolation;

            h = computationalGrid.h;
            [nL,nR,nB,nT,widths] = pml.cellWidths(h);
            nPhysical = computationalGrid.N;
            nExtended = nPhysical+[nL+nR,nB+nT];
            extendedDomain = [ ...
                computationalGrid.dom(1,1)-widths(1), ...
                computationalGrid.dom(1,2)+widths(2); ...
                computationalGrid.dom(2,1)-widths(3), ...
                computationalGrid.dom(2,2)+widths(4)];

            xExtended = Grid1DUnifPrimal(extendedDomain(1,:),nExtended(1));
            yExtended = Grid1DUnifPrimal(extendedDomain(2,:),nExtended(2));
            obj.extendedGrid = GridSpace2D(xExtended,yExtended);
            obj.extendedVelocityXGrid = GridSpace2D( ...
                Grid1DUnifDual(xExtended),yExtended);
            obj.extendedVelocityYGrid = GridSpace2D( ...
                xExtended,Grid1DUnifDual(yExtended));

            obj.physicalPressureXIndices = nL+(1:nPhysical(1));
            obj.physicalPressureYIndices = nB+(1:nPhysical(2));
            obj.physicalVelocityXXIndices = nL+(1:nPhysical(1)-1);
            obj.physicalVelocityXYIndices = nB+(1:nPhysical(2));
            obj.physicalVelocityYXIndices = nL+(1:nPhysical(1));
            obj.physicalVelocityYYIndices = nB+(1:nPhysical(2)-1);

            [Kphysical,BxPhysical,ByPhysical,~,obj.cMax] = ...
                sampleMedium2D(computationalGrid,medium);
            obj.K = GridFunction(obj.extendedGrid, ...
                Label=Kphysical.label,Units=Kphysical.units);
            obj.K.setValues(obj.extendConstant2D( ...
                Kphysical.values,nL,nR,nB,nT));
            obj.Bx = GridFunction(obj.extendedVelocityXGrid, ...
                Label=BxPhysical.label,Units=BxPhysical.units);
            obj.Bx.setValues(obj.extendConstant2D( ...
                BxPhysical.values,nL,nR,nB,nT));
            obj.By = GridFunction(obj.extendedVelocityYGrid, ...
                Label=ByPhysical.label,Units=ByPhysical.units);
            obj.By.setValues(obj.extendConstant2D( ...
                ByPhysical.values,nL,nR,nB,nT));

            sigmaXP = pml.dampingX(obj.extendedGrid.x.pts, ...
                computationalGrid.dom(1,:),widths(1:2),obj.cMax);
            sigmaYP = pml.dampingY(obj.extendedGrid.y.pts, ...
                computationalGrid.dom(2,:),widths(3:4),obj.cMax);
            obj.sigmaPX = repmat(sigmaXP(:),1,obj.extendedGrid.N(2));
            obj.sigmaPY = repmat(sigmaYP(:).',obj.extendedGrid.N(1),1);

            sigmaXV = pml.dampingX(obj.extendedVelocityXGrid.x.pts, ...
                computationalGrid.dom(1,:),widths(1:2),obj.cMax);
            sigmaYV = pml.dampingY(obj.extendedVelocityYGrid.y.pts, ...
                computationalGrid.dom(2,:),widths(3:4),obj.cMax);
            obj.sigmaVX = repmat( ...
                sigmaXV(:),1,obj.extendedVelocityXGrid.N(2));
            obj.sigmaVY = repmat( ...
                sigmaYV(:).',obj.extendedVelocityYGrid.N(1),1);

            obj.cfl = obj.cMax*obj.dt*sqrt(sum(1./h.^2));
            cflLimit = obj.stabilityLimit();
            if options.EnforceCFL && obj.cfl > cflLimit+100*eps
                error('AcousticSolver2D_PML:CFLViolation', ...
                    ['The multidimensional CFL number %.6g exceeds the ' ...
                     'order-%d stability limit %.6g.'], ...
                    obj.cfl,obj.spatialOrder,cflLimit);
            end

            if obj.velocityExtrapolation == "none"
                obj.validateVelocityOutputDomain( ...
                    obj.extendedVelocityXGrid,"x-velocity");
                obj.validateVelocityOutputDomain( ...
                    obj.extendedVelocityYGrid,"y-velocity");
            end
        end

        function [pressure,velocityX,velocityY] = solve(obj,p0,vx0,vy0,options)
            arguments
                obj
                p0
                vx0
                vy0
                options.Source (1,1) AcousticSource2D = ...
                    AcousticSource2D.zero()
            end

            pressure = GridFunctionTimeSeries(obj.outputGrid,obj.outputTimes, ...
                Label="Pressure",Units="Pa",TimeUnits="s");
            velocityX = GridFunctionTimeSeries(obj.outputGrid,obj.outputTimes, ...
                Label="x-particle velocity",Units="m/s",TimeUnits="s");
            velocityY = GridFunctionTimeSeries(obj.outputGrid,obj.outputTimes, ...
                Label="y-particle velocity",Units="m/s",TimeUnits="s");

            dt = obj.dt;
            dx = obj.computationalGrid.h(1);
            dy = obj.computationalGrid.h(2);
            times = obj.computationalTimes;
            pressureGrid = obj.extendedGrid;
            velocityXGrid = obj.extendedVelocityXGrid;
            velocityYGrid = obj.extendedVelocityYGrid;

            pInitial = obj.initialDataOnExtendedGrid( ...
                p0,pressureGrid,obj.computationalGrid, ...
                obj.physicalPressureXIndices, ...
                obj.physicalPressureYIndices,'pressure');
            px = 0.5*pInitial;
            py = 0.5*pInitial;
            [px,py] = obj.applyPressureBoundary(px,py);

            physicalVXGrid = GridSpace2D( ...
                Grid1DUnifDual(obj.computationalGrid.x), ...
                obj.computationalGrid.y);
            physicalVYGrid = GridSpace2D( ...
                obj.computationalGrid.x, ...
                Grid1DUnifDual(obj.computationalGrid.y));
            vx0Values = obj.initialDataOnExtendedGrid( ...
                vx0,velocityXGrid,physicalVXGrid, ...
                obj.physicalVelocityXXIndices, ...
                obj.physicalVelocityXYIndices,'x-velocity');
            vy0Values = obj.initialDataOnExtendedGrid( ...
                vy0,velocityYGrid,physicalVYGrid, ...
                obj.physicalVelocityYXIndices, ...
                obj.physicalVelocityYYIndices,'y-velocity');

            pressureSource = options.Source.preparePressure( ...
                obj.computationalGrid);
            velocityXSource = options.Source.prepareVelocityX(physicalVXGrid);
            velocityYSource = options.Source.prepareVelocityY(physicalVYGrid);

            [aPX,bPX] = obj.dampingCoefficients(obj.sigmaPX,dt);
            [aPY,bPY] = obj.dampingCoefficients(obj.sigmaPY,dt);
            [aVX,bVX] = obj.dampingCoefficients(obj.sigmaVX,dt);
            [aVY,bVY] = obj.dampingCoefficients(obj.sigmaVY,dt);
            [aVXHalf,bVXHalf] = obj.dampingCoefficients(obj.sigmaVX,dt/2);
            [aVYHalf,bVYHalf] = obj.dampingCoefficients(obj.sigmaVY,dt/2);

            pOutputIndex = 1;
            vxOutputIndex = 1;
            vyOutputIndex = 1;
            if obj.timesEqual(obj.outputTimes(1),times(1))
                p = px+py;
                pressure.setSnapshot(1,obj.interpolateSpace( ...
                    p,pressureGrid,false));
                velocityX.setSnapshot(1,obj.interpolateSpace( ...
                    vx0Values,velocityXGrid,true));
                velocityY.setSnapshot(1,obj.interpolateSpace( ...
                    vy0Values,velocityYGrid,true));
                pOutputIndex = 2;
                vxOutputIndex = 2;
                vyOutputIndex = 2;
            end

            p = px+py;
            DxP = obj.derivativePrimalToDualX(p,dx);
            DyP = obj.derivativePrimalToDualY(p,dy);
            fx = obj.embedVelocityXSource(options.Source.evaluatePrepared( ...
                velocityXSource,times(1)));
            fy = obj.embedVelocityYSource(options.Source.evaluatePrepared( ...
                velocityYSource,times(1)));
            vxHalf = aVXHalf.*vx0Values ...
                -bVXHalf.*obj.Bx.values.*DxP+bVXHalf.*fx;
            vyHalf = aVYHalf.*vy0Values ...
                -bVYHalf.*obj.By.values.*DyP+bVYHalf.*fy;
            tVRight = times(1)+dt/2;

            vxOutputIndex = obj.storeBracketedOutputs( ...
                velocityX,vxOutputIndex,vx0Values,times(1), ...
                vxHalf,tVRight,velocityXGrid,true);
            vyOutputIndex = obj.storeBracketedOutputs( ...
                velocityY,vyOutputIndex,vy0Values,times(1), ...
                vyHalf,tVRight,velocityYGrid,true);

            for n = 2:obj.Nt
                pOld = px+py;
                DxVx = obj.derivativeDualToPrimalX(vxHalf,dx);
                DyVy = obj.derivativeDualToPrimalY(vyHalf,dy);
                fp = obj.embedPressureSource( ...
                    options.Source.evaluatePrepared(pressureSource, ...
                    0.5*(times(n-1)+times(n))));
                interiorX = 2:size(px,1)-1;
                interiorY = 2:size(px,2)-1;

                px(interiorX,interiorY) = ...
                    aPX(interiorX,interiorY).*px(interiorX,interiorY) ...
                    -bPX(interiorX,interiorY).* ...
                    obj.K.values(interiorX,interiorY).* ...
                    DxVx(:,interiorY) ...
                    +0.5*bPX(interiorX,interiorY).* ...
                    fp(interiorX,interiorY);
                py(interiorX,interiorY) = ...
                    aPY(interiorX,interiorY).*py(interiorX,interiorY) ...
                    -bPY(interiorX,interiorY).* ...
                    obj.K.values(interiorX,interiorY).* ...
                    DyVy(interiorX,:) ...
                    +0.5*bPY(interiorX,interiorY).* ...
                    fp(interiorX,interiorY);
                [px,py] = obj.applyPressureBoundary(px,py);
                p = px+py;

                pOutputIndex = obj.storeBracketedOutputs( ...
                    pressure,pOutputIndex,pOld,times(n-1),p,times(n), ...
                    pressureGrid,false);

                vxOld = vxHalf;
                vyOld = vyHalf;
                oldVelocityTime = tVRight;
                DxP = obj.derivativePrimalToDualX(p,dx);
                DyP = obj.derivativePrimalToDualY(p,dy);
                fx = obj.embedVelocityXSource( ...
                    options.Source.evaluatePrepared(velocityXSource,times(n)));
                fy = obj.embedVelocityYSource( ...
                    options.Source.evaluatePrepared(velocityYSource,times(n)));
                vxHalf = aVX.*vxOld ...
                    -bVX.*obj.Bx.values.*DxP+bVX.*fx;
                vyHalf = aVY.*vyOld ...
                    -bVY.*obj.By.values.*DyP+bVY.*fy;
                tVRight = times(n)+dt/2;

                vxOutputIndex = obj.storeBracketedOutputs( ...
                    velocityX,vxOutputIndex,vxOld,oldVelocityTime, ...
                    vxHalf,tVRight,velocityXGrid,true);
                vyOutputIndex = obj.storeBracketedOutputs( ...
                    velocityY,vyOutputIndex,vyOld,oldVelocityTime, ...
                    vyHalf,tVRight,velocityYGrid,true);
            end

            if pOutputIndex <= obj.NtOutput || ...
                    vxOutputIndex <= obj.NtOutput || ...
                    vyOutputIndex <= obj.NtOutput
                error('AcousticSolver2D_PML:UnfilledOutput', ...
                    'One or more requested output times were not generated.');
            end
        end
    end

    methods (Access = private)
        function limit = stabilityLimit(obj)
            if obj.spatialOrder == 2
                limit = 1;
            else
                limit = 6/7;
            end
        end

        function [a,b] = dampingCoefficients(~,sigma,timeStep)
            a = (1-0.5*timeStep*sigma)./(1+0.5*timeStep*sigma);
            b = timeStep./(1+0.5*timeStep*sigma);
        end

        function DxP = derivativePrimalToDualX(obj,p,dx)
            if obj.spatialOrder == 2
                DxP = diff(p,1,1)/dx;
                return
            end
            pExtended = [-p(2,:);p;-p(end-1,:)];
            i = (1:size(p,1)-1).';
            DxP = (pExtended(i,:)-27*pExtended(i+1,:) ...
                +27*pExtended(i+2,:)-pExtended(i+3,:))/(24*dx);
        end

        function DyP = derivativePrimalToDualY(obj,p,dy)
            if obj.spatialOrder == 2
                DyP = diff(p,1,2)/dy;
                return
            end
            pExtended = [-p(:,2),p,-p(:,end-1)];
            j = 1:size(p,2)-1;
            DyP = (pExtended(:,j)-27*pExtended(:,j+1) ...
                +27*pExtended(:,j+2)-pExtended(:,j+3))/(24*dy);
        end

        function DxV = derivativeDualToPrimalX(obj,v,dx)
            if obj.spatialOrder == 2
                DxV = diff(v,1,1)/dx;
                return
            end
            vExtended = [v(1,:);v;v(end,:)];
            i = (1:size(v,1)-1).';
            DxV = (vExtended(i,:)-27*vExtended(i+1,:) ...
                +27*vExtended(i+2,:)-vExtended(i+3,:))/(24*dx);
        end

        function DyV = derivativeDualToPrimalY(obj,v,dy)
            if obj.spatialOrder == 2
                DyV = diff(v,1,2)/dy;
                return
            end
            vExtended = [v(:,1),v,v(:,end)];
            j = 1:size(v,2)-1;
            DyV = (vExtended(:,j)-27*vExtended(:,j+1) ...
                +27*vExtended(:,j+2)-vExtended(:,j+3))/(24*dy);
        end

        function values = initialDataOnExtendedGrid(obj,data, ...
                extendedGrid,physicalGrid,physicalX,physicalY,fieldName)
            if isa(data,'function_handle')
                target = GridFunction(extendedGrid);
                try
                    target.setValues(data);
                catch cause
                    error('AcousticSolver2D_PML:InvalidInitialData', ...
                        'Invalid initial %s data: %s',fieldName,cause.message);
                end
                values = target.values;
                return
            end

            physicalValues = obj.dataOnGrid(data,physicalGrid,fieldName);
            nL = physicalX(1)-1;
            nR = extendedGrid.N(1)-physicalX(end);
            nB = physicalY(1)-1;
            nT = extendedGrid.N(2)-physicalY(end);
            values = obj.extendConstant2D(physicalValues,nL,nR,nB,nT);
        end

        function values = dataOnGrid(~,data,grid,fieldName)
            target = GridFunction(grid);
            if isa(data,'GridFunction')
                if data.dim ~= 2
                    error('AcousticSolver2D_PML:DimensionMismatch', ...
                        'Initial %s data must be two-dimensional.',fieldName);
                end
                if isequal(data.grid,grid)
                    values = data.values;
                else
                    target.interpolateFrom(data);
                    values = target.values;
                end
            else
                try
                    target.setValues(data);
                    values = target.values;
                catch cause
                    error('AcousticSolver2D_PML:InvalidInitialData', ...
                        'Invalid initial %s data: %s',fieldName,cause.message);
                end
            end
        end

        function values = extendConstant2D(~,values,nL,nR,nB,nT)
            ix = [ones(1,nL),1:size(values,1), ...
                size(values,1)*ones(1,nR)];
            iy = [ones(1,nB),1:size(values,2), ...
                size(values,2)*ones(1,nT)];
            values = values(ix,iy);
        end

        function values = embedPressureSource(obj,physicalValues)
            values = zeros(obj.extendedGrid.N);
            values(obj.physicalPressureXIndices, ...
                obj.physicalPressureYIndices) = physicalValues;
        end

        function values = embedVelocityXSource(obj,physicalValues)
            values = zeros(obj.extendedVelocityXGrid.N);
            values(obj.physicalVelocityXXIndices, ...
                obj.physicalVelocityXYIndices) = physicalValues;
        end

        function values = embedVelocityYSource(obj,physicalValues)
            values = zeros(obj.extendedVelocityYGrid.N);
            values(obj.physicalVelocityYXIndices, ...
                obj.physicalVelocityYYIndices) = physicalValues;
        end

        function [px,py] = applyPressureBoundary(~,px,py)
            px([1,end],:) = 0;
            px(:,[1,end]) = 0;
            py([1,end],:) = 0;
            py(:,[1,end]) = 0;
        end

        function nextIndex = storeBracketedOutputs(obj,series,nextIndex, ...
                leftValues,leftTime,rightValues,rightTime,sourceGrid,isVelocity)
            while nextIndex <= obj.NtOutput && ...
                    obj.outputTimes(nextIndex) <= rightTime+obj.timeTolerance()
                outputTime = obj.outputTimes(nextIndex);
                if outputTime < leftTime-obj.timeTolerance()
                    error('AcousticSolver2D_PML:OutputOrdering', ...
                        'Output time was not bracketed by computational states.');
                end
                theta = (outputTime-leftTime)/(rightTime-leftTime);
                theta = min(max(theta,0),1);
                values = (1-theta)*leftValues+theta*rightValues;
                series.setSnapshot(nextIndex,obj.interpolateSpace( ...
                    values,sourceGrid,isVelocity));
                nextIndex = nextIndex+1;
            end
        end

        function values = interpolateSpace(obj,data,sourceGrid,isVelocity)
            if isequal(sourceGrid,obj.outputGrid)
                values = data;
            else
                if isVelocity && obj.velocityExtrapolation == "linear"
                    extrapolation = "linear";
                else
                    extrapolation = "none";
                end
                interpolant = griddedInterpolant( ...
                    {sourceGrid.x.pts,sourceGrid.y.pts},data, ...
                    "linear",extrapolation);
                [X,Y] = obj.outputGrid.mesh();
                values = interpolant(X,Y);
            end
            if any(~isfinite(values),'all')
                error('AcousticSolver2D_PML:SpatialInterpolationFailed', ...
                    'Spatial interpolation produced a nonfinite value.');
            end
        end

        function validateVelocityOutputDomain(obj,velocityGrid,component)
            if any(obj.outputGrid.dom(:,1) < velocityGrid.dom(:,1)) || ...
                    any(obj.outputGrid.dom(:,2) > velocityGrid.dom(:,2))
                error('AcousticSolver2D_PML:VelocityExtrapolationRequired', ...
                    ['The output grid extends beyond the extended %s grid. ' ...
                     'Use VelocityExtrapolation="linear" or select interior ' ...
                     'output points.'],component);
            end
        end

        function value = timesEqual(obj,a,b)
            value = abs(a-b) <= obj.timeTolerance();
        end

        function value = timeTolerance(obj)
            value = 100*eps(max(1,max(abs(obj.computationalTimes))));
        end
    end
end
