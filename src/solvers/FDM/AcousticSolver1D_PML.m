classdef AcousticSolver1D_PML < AcousticSolver
    % AcousticSolver1D_PML solves 1D acoustics with two PML layers.
    % Pressure uses an extended uniform primal grid and velocity its dual.
    % SpatialOrder may be 2 or 4; leapfrog time stepping is second order.

    properties (SetAccess = protected)
        medium
        pml
        spatialOrder
        extendedGrid
        extendedDualGrid
        physicalPressureIndices
        physicalVelocityIndices
        K
        B
        cMax
        cfl
        sigmaP
        sigmaV
        velocityExtrapolation
    end

    methods
        function obj = AcousticSolver1D_PML( ...
                computationalGrid,computationalTimes, ...
                outputGrid,outputTimes,medium,pml,options)
            arguments
                computationalGrid (1,1) GridSpace1D
                computationalTimes double {mustBeFinite}
                outputGrid (1,1) GridSpace1D
                outputTimes double {mustBeFinite}
                medium (1,1) AcousticMedium1D
                pml (1,1) PML1D
                options.SpatialOrder (1,1) double ...
                    {mustBeMember(options.SpatialOrder,[2,4])} = 2
                options.EnforceCFL (1,1) logical = true
                options.VelocityExtrapolation (1,1) string ...
                    {mustBeMember(options.VelocityExtrapolation, ...
                    ["linear","none"])} = "linear"
            end

            if ~isa(computationalGrid.x,'Grid1DUnifPrimal')
                error('AcousticSolver1D_PML:PrimalGridRequired', ...
                    'The physical pressure grid must be uniform and primal.');
            end
            if outputGrid.dom(1) < computationalGrid.dom(1) || ...
                    outputGrid.dom(2) > computationalGrid.dom(2)
                error('AcousticSolver1D_PML:OutputGridOutsideDomain', ...
                    'Output points must lie in the physical domain.');
            end

            obj@AcousticSolver(computationalGrid,computationalTimes, ...
                outputGrid,outputTimes);

            obj.medium = medium;
            obj.pml = pml;
            obj.spatialOrder = options.SpatialOrder;
            obj.velocityExtrapolation = options.VelocityExtrapolation;

            dx = computationalGrid.h;
            [nLeft,nRight,widths] = pml.cellWidths(dx);
            nPhysical = computationalGrid.x.N;
            nExtended = nPhysical+nLeft+nRight;
            extendedDomain = [computationalGrid.dom(1)-widths(1), ...
                computationalGrid.dom(2)+widths(2)];

            obj.extendedGrid = GridSpace1D( ...
                Grid1DUnifPrimal(extendedDomain,nExtended));
            obj.extendedDualGrid = GridSpace1D( ...
                Grid1DUnifDual(obj.extendedGrid.x));
            obj.physicalPressureIndices = nLeft+(1:nPhysical);
            obj.physicalVelocityIndices = nLeft+(1:nPhysical-1);

            [Kphysical,Bphysical,~,obj.cMax] = sampleMedium1D( ...
                computationalGrid,medium,BetaGrid="dual");
            obj.K = GridFunction(obj.extendedGrid, ...
                Label=Kphysical.label,Units=Kphysical.units);
            obj.K.setValues(obj.extendConstant( ...
                Kphysical.values,nLeft,nRight));
            obj.B = GridFunction(obj.extendedDualGrid, ...
                Label=Bphysical.label,Units=Bphysical.units);
            obj.B.setValues(obj.extendConstant( ...
                Bphysical.values,nLeft,nRight));

            obj.sigmaP = pml.damping(obj.extendedGrid.x.pts, ...
                computationalGrid.dom,widths,obj.cMax);
            obj.sigmaV = pml.damping(obj.extendedDualGrid.x.pts, ...
                computationalGrid.dom,widths,obj.cMax);

            obj.cfl = obj.cMax*obj.dt/dx;
            cflLimit = 1;
            if obj.spatialOrder == 4
                cflLimit = 6/7;
            end
            if options.EnforceCFL && obj.cfl > cflLimit+100*eps
                error('AcousticSolver1D_PML:CFLViolation', ...
                    ['The CFL number %.6g exceeds the order-%d ' ...
                     'stability limit %.6g.'], ...
                    obj.cfl,obj.spatialOrder,cflLimit);
            end

            dualPoints = obj.extendedDualGrid.x.pts;
            outputPoints = outputGrid.x.pts;
            if obj.velocityExtrapolation == "none" && ...
                    (outputPoints(1) < dualPoints(1) || ...
                     outputPoints(end) > dualPoints(end))
                error('AcousticSolver1D_PML:VelocityExtrapolationRequired', ...
                    ['Use VelocityExtrapolation="linear" or choose ' ...
                     'output points inside the extended velocity grid.']);
            end
        end

        function [pressure,velocity] = solve(obj,p0,v0,options)
            arguments
                obj
                p0
                v0
                options.Source (1,1) AcousticSource1D = ...
                    AcousticSource1D.zero()
            end

            pressure = GridFunctionTimeSeries( ...
                obj.outputGrid,obj.outputTimes, ...
                Label="Pressure",Units="Pa",TimeUnits="s");
            velocity = GridFunctionTimeSeries( ...
                obj.outputGrid,obj.outputTimes, ...
                Label="Particle velocity",Units="m/s",TimeUnits="s");

            dt = obj.dt;
            dx = obj.computationalGrid.h;
            times = obj.computationalTimes;
            primalGrid = obj.extendedGrid;
            dualGrid = obj.extendedDualGrid;

            p = obj.initialDataOnExtendedGrid( ...
                p0,primalGrid,obj.computationalGrid, ...
                obj.physicalPressureIndices,'pressure');
            v0Values = obj.initialDataOnExtendedGrid( ...
                v0,dualGrid,GridSpace1D( ...
                Grid1DUnifDual(obj.computationalGrid.x)), ...
                obj.physicalVelocityIndices,'velocity');

            % Only the outer endpoints are boundaries. The original-domain
            % endpoints are interior points at the PML interfaces.
            p([1,end]) = 0;

            % Sources are constructed only on the physical grids. Their
            % values are embedded into the extended arrays below, so no
            % source term can be injected into either PML layer.
            physicalPressureGrid = obj.computationalGrid;
            physicalVelocityGrid = GridSpace1D( ...
                Grid1DUnifDual(obj.computationalGrid.x));
            pressureSource = options.Source.preparePressure( ...
                physicalPressureGrid);
            velocitySource = options.Source.prepareVelocity( ...
                physicalVelocityGrid);

            aP = (1-0.5*dt*obj.sigmaP)./(1+0.5*dt*obj.sigmaP);
            bP = dt./(1+0.5*dt*obj.sigmaP);
            aV = (1-0.5*dt*obj.sigmaV)./(1+0.5*dt*obj.sigmaV);
            bV = dt./(1+0.5*dt*obj.sigmaV);

            pOutputIndex = 1;
            vOutputIndex = 1;
            if obj.timesEqual(obj.outputTimes(1),times(1))
                pressure.setSnapshot(1,obj.interpolateSpace( ...
                    p,primalGrid,false));
                velocity.setSnapshot(1,obj.interpolateSpace( ...
                    v0Values,dualGrid,true));
                pOutputIndex = 2;
                vOutputIndex = 2;
            end

            % Trapezoidal damping over the initial half step.
            aVHalf = (1-0.25*dt*obj.sigmaV) ...
                ./(1+0.25*dt*obj.sigmaV);
            bVHalf = (0.5*dt)./(1+0.25*dt*obj.sigmaV);
            DxP = obj.derivativePrimalToDual(p,dx);
            velocityForcing = zeros(dualGrid.N,1);
            velocityForcing(obj.physicalVelocityIndices) = ...
                options.Source.evaluatePrepared(velocitySource,times(1));
            vHalf = aVHalf.*v0Values ...
                - bVHalf.*obj.B.values.*DxP ...
                + bVHalf.*velocityForcing;
            tVRight = times(1)+dt/2;
            vOutputIndex = obj.storeBracketedOutputs( ...
                velocity,vOutputIndex,v0Values,times(1), ...
                vHalf,tVRight,dualGrid,true);

            for n = 2:obj.Nt
                pOld = p;
                DxV = obj.derivativeDualToPrimal(vHalf,dx);
                tHalf = 0.5*(times(n-1)+times(n));
                pressureForcing = zeros(primalGrid.N,1);
                pressureForcing(obj.physicalPressureIndices) = ...
                    options.Source.evaluatePrepared(pressureSource,tHalf);
                interior = 2:numel(p)-1;
                p(interior) = aP(interior).*pOld(interior) ...
                    - bP(interior).*obj.K.values(interior).*DxV ...
                    + bP(interior).*pressureForcing(interior);
                p([1,end]) = 0;

                pOutputIndex = obj.storeBracketedOutputs( ...
                    pressure,pOutputIndex,pOld,times(n-1),p,times(n), ...
                    primalGrid,false);

                vOld = vHalf;
                DxP = obj.derivativePrimalToDual(p,dx);
                velocityForcing = zeros(dualGrid.N,1);
                velocityForcing(obj.physicalVelocityIndices) = ...
                    options.Source.evaluatePrepared(velocitySource,times(n));
                vHalf = aV.*vOld ...
                    - bV.*obj.B.values.*DxP ...
                    + bV.*velocityForcing;
                oldVelocityTime = tVRight;
                tVRight = times(n)+dt/2;
                vOutputIndex = obj.storeBracketedOutputs( ...
                    velocity,vOutputIndex,vOld,oldVelocityTime, ...
                    vHalf,tVRight,dualGrid,true);
            end

            if pOutputIndex <= obj.NtOutput || vOutputIndex <= obj.NtOutput
                error('AcousticSolver1D_PML:UnfilledOutput', ...
                    'One or more requested output times were not generated.');
            end
        end
    end

    methods (Access = private)
        function DxP = derivativePrimalToDual(obj,p,dx)
            if obj.spatialOrder == 2
                DxP = diff(p)/dx;
                return
            end
            nPressure = numel(p);
            pExtended = [-p(2);p;-p(end-1)];
            i = (1:nPressure-1).';
            DxP = (pExtended(i)-27*pExtended(i+1) ...
                +27*pExtended(i+2)-pExtended(i+3))/(24*dx);
        end

        function DxV = derivativeDualToPrimal(obj,v,dx)
            if obj.spatialOrder == 2
                DxV = diff(v)/dx;
                return
            end
            nVelocity = numel(v);
            vExtended = [v(1);v;v(end)];
            i = (1:nVelocity-1).';
            DxV = (vExtended(i)-27*vExtended(i+1) ...
                +27*vExtended(i+2)-vExtended(i+3))/(24*dx);
        end

        function values = initialDataOnExtendedGrid(obj,data, ...
                extendedGrid,physicalGrid,physicalIndices,fieldName)
            if isa(data,'function_handle')
                target = GridFunction(extendedGrid);
                try
                    target.setValues(data);
                catch cause
                    error('AcousticSolver1D_PML:InvalidInitialData', ...
                        'Invalid initial %s data: %s',fieldName,cause.message);
                end
                values = target.values;
                return
            end

            physicalValues = obj.dataOnGrid(data,physicalGrid,fieldName);
            nLeft = physicalIndices(1)-1;
            nRight = extendedGrid.N-physicalIndices(end);
            values = obj.extendConstant(physicalValues,nLeft,nRight);
        end

        function values = dataOnGrid(~,data,grid,fieldName)
            target = GridFunction(grid);
            if isa(data,'GridFunction')
                if data.dim ~= grid.dim
                    error('AcousticSolver1D_PML:DimensionMismatch', ...
                        'Initial %s data must be one-dimensional.',fieldName);
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
                    error('AcousticSolver1D_PML:InvalidInitialData', ...
                        'Invalid initial %s data: %s',fieldName,cause.message);
                end
            end
        end

        function values = extendConstant(~,values,nLeft,nRight)
            values = values(:);
            values = [repmat(values(1),nLeft,1);values; ...
                repmat(values(end),nRight,1)];
        end

        function nextIndex = storeBracketedOutputs(obj,series,nextIndex, ...
                leftValues,leftTime,rightValues,rightTime,sourceGrid,isVelocity)
            while nextIndex <= obj.NtOutput && ...
                    obj.outputTimes(nextIndex) <= rightTime+obj.timeTolerance()
                outputTime = obj.outputTimes(nextIndex);
                if outputTime < leftTime-obj.timeTolerance()
                    error('AcousticSolver1D_PML:OutputOrdering', ...
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
            sourcePoints = sourceGrid.x.pts;
            outputPoints = obj.outputGrid.x.pts;
            if isequal(sourceGrid,obj.outputGrid)
                values = data;
            elseif isVelocity && obj.velocityExtrapolation == "linear"
                values = interp1(sourcePoints,data,outputPoints,'linear','extrap');
            else
                values = interp1(sourcePoints,data,outputPoints,'linear');
            end
            if any(~isfinite(values),'all')
                error('AcousticSolver1D_PML:SpatialInterpolationFailed', ...
                    'Spatial interpolation produced a nonfinite value.');
            end
            values = values(:);
        end

        function value = timesEqual(obj,a,b)
            value = abs(a-b) <= obj.timeTolerance();
        end

        function value = timeTolerance(obj)
            value = 100*eps(max(1,max(abs(obj.computationalTimes))));
        end
    end
end
