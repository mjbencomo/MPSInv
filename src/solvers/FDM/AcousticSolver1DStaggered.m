classdef AcousticSolver1DStaggered < AcousticSolver
    % AcousticSolver1DStaggered solves the 1D first-order acoustic system.
    % Pressure is stored on a uniform primal grid and particle velocity on
    % its dual grid. SpatialOrder selects second- or fourth-order centered
    % staggered differences. Leapfrog time stepping is second order.

    properties (SetAccess = protected)
        medium
        K
        B
        cMax
        cfl
        spatialOrder
        velocityExtrapolation
    end

    methods
        function obj = AcousticSolver1DStaggered( ...
                computationalGrid,computationalTimes, ...
                outputGrid,outputTimes,medium,options)
            arguments
                computationalGrid (1,1) GridSpace1D
                computationalTimes double {mustBeFinite}
                outputGrid (1,1) GridSpace1D
                outputTimes double {mustBeFinite}
                medium (1,1) AcousticMedium1D
                options.SpatialOrder (1,1) double ...
                    {mustBeMember(options.SpatialOrder,[2,4])} = 2
                options.EnforceCFL (1,1) logical = true
                options.VelocityExtrapolation (1,1) string ...
                    {mustBeMember(options.VelocityExtrapolation, ...
                    ["linear","none"])} = "linear"
            end

            if ~isa(computationalGrid.x,'Grid1DUnifPrimal')
                error('AcousticSolver1DStaggered:PrimalGridRequired', ...
                    'Computational pressure requires a uniform primal grid.');
            end

            if options.SpatialOrder == 4 && computationalGrid.x.N < 5
                error('AcousticSolver1DStaggered:GridTooSmall', ...
                    ['At least five pressure grid points are required for ' ...
                     'the fourth-order staggered stencil.']);
            end

            if outputGrid.dom(1) < computationalGrid.dom(1) || ...
                    outputGrid.dom(2) > computationalGrid.dom(2)
                error('AcousticSolver1DStaggered:OutputGridOutsideDomain', ...
                    ['Output points must lie within the computational ' ...
                     'spatial domain.']);
            end

            obj@AcousticSolver(computationalGrid,computationalTimes, ...
                outputGrid,outputTimes);

            obj.medium = medium;
            obj.spatialOrder = options.SpatialOrder;
            [obj.K,obj.B,~,obj.cMax] = sampleMedium1D( ...
                computationalGrid,medium,BetaGrid="dual");
            obj.cfl = obj.cMax*obj.dt/computationalGrid.h;
            obj.velocityExtrapolation = options.VelocityExtrapolation;

            cflLimit = obj.stabilityLimit();
            if options.EnforceCFL && obj.cfl > cflLimit+100*eps
                error('AcousticSolver1DStaggered:CFLViolation', ...
                    ['The CFL number cMax*dt/dx = %.6g exceeds the ' ...
                     'order-%d stability limit %.6g.'], ...
                    obj.cfl,obj.spatialOrder,cflLimit);
            end

            dualPoints = obj.B.grid.x.pts;
            outputPoints = outputGrid.x.pts;
            if obj.velocityExtrapolation == "none" && ...
                    (outputPoints(1) < dualPoints(1) || ...
                     outputPoints(end) > dualPoints(end))
                error('AcousticSolver1DStaggered:VelocityExtrapolationRequired', ...
                    ['The output grid extends beyond the computational ' ...
                     'velocity grid. Use VelocityExtrapolation="linear" ' ...
                     'or choose interior output points.']);
            end
        end

        function [pressure,velocity] = solve(obj,p0,v0,options)
            % Solve and sample both fields on the requested output grid.
            % Homogeneous Dirichlet conditions are imposed on pressure.
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

            primalGrid = obj.computationalGrid;
            dualGrid = obj.B.grid;
            dx = primalGrid.h;
            dt = obj.dt;
            times = obj.computationalTimes;

            % Sample all spatial source terms once. During time stepping,
            % only their scalar time factors and a matrix-vector product
            % are evaluated. All-localized sources retain sparse storage.
            pressureSource = options.Source.preparePressure(primalGrid);
            velocitySource = options.Source.prepareVelocity(dualGrid);

            p = obj.dataOnGrid(p0,primalGrid,'pressure');
            v0Values = obj.dataOnGrid(v0,dualGrid,'velocity');
            p([1,end]) = 0;

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

            % Leapfrog start: v(t_0) -> v(t_{1/2}).
            DxP = obj.derivativePrimalToDual(p,dx);
            velocityForcing = options.Source.evaluatePrepared( ...
                velocitySource,times(1));
            vHalf = v0Values ...
                -(dt/2)*obj.B.values.*DxP ...
                +(dt/2)*velocityForcing;
            tVRight = times(1)+dt/2;
            vOutputIndex = obj.storeBracketedOutputs( ...
                velocity,vOutputIndex,v0Values,times(1), ...
                vHalf,tVRight,dualGrid,true);

            for n = 2:obj.Nt
                pOld = p;
                DxV = obj.derivativeDualToPrimal(vHalf,dx);
                tHalf = 0.5*(times(n-1)+times(n));
                pressureForcing = options.Source.evaluatePrepared( ...
                    pressureSource,tHalf);

                p(2:end-1) = pOld(2:end-1) ...
                    -dt*obj.K.values(2:end-1).*DxV ...
                    +dt*pressureForcing(2:end-1);
                p([1,end]) = 0;

                pOutputIndex = obj.storeBracketedOutputs( ...
                    pressure,pOutputIndex,pOld,times(n-1),p,times(n), ...
                    primalGrid,false);

                vOld = vHalf;
                oldVelocityTime = tVRight;
                DxP = obj.derivativePrimalToDual(p,dx);
                velocityForcing = options.Source.evaluatePrepared( ...
                    velocitySource,times(n));
                vHalf = vOld ...
                    -dt*obj.B.values.*DxP ...
                    +dt*velocityForcing;
                tVRight = times(n)+dt/2;

                vOutputIndex = obj.storeBracketedOutputs( ...
                    velocity,vOutputIndex,vOld,oldVelocityTime, ...
                    vHalf,tVRight,dualGrid,true);
            end

            if pOutputIndex <= obj.NtOutput || vOutputIndex <= obj.NtOutput
                error('AcousticSolver1DStaggered:UnfilledOutput', ...
                    'One or more requested output times were not generated.');
            end
        end
    end

    methods (Access = private)
        function limit = stabilityLimit(obj)
            if obj.spatialOrder == 2
                limit = 1;
            else
                % Maximum fourth-order staggered derivative symbol is
                % 7/(3*dx), which gives CFL <= 6/7 for leapfrog.
                limit = 6/7;
            end
        end

        function DxP = derivativePrimalToDual(obj,p,dx)
            if obj.spatialOrder == 2
                DxP = diff(p)/dx;
                return
            end

            % Odd pressure extension preserves p=0 at both boundaries.
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

            % Even velocity extension supplies the fourth-order ghosts.
            nVelocity = numel(v);
            vExtended = [v(1);v;v(end)];
            i = (1:nVelocity-1).';
            DxV = (vExtended(i)-27*vExtended(i+1) ...
                +27*vExtended(i+2)-vExtended(i+3))/(24*dx);
        end

        function nextIndex = storeBracketedOutputs(obj,series,nextIndex, ...
                leftValues,leftTime,rightValues,rightTime,sourceGrid,isVelocity)
            while nextIndex <= obj.NtOutput && ...
                    obj.outputTimes(nextIndex) <= rightTime+obj.timeTolerance()
                outputTime = obj.outputTimes(nextIndex);
                if outputTime < leftTime-obj.timeTolerance()
                    error('AcousticSolver1DStaggered:OutputOrdering', ...
                        'Output time was not bracketed by computational states.');
                end

                theta = (outputTime-leftTime)/(rightTime-leftTime);
                theta = min(max(theta,0),1);
                values = (1-theta)*leftValues+theta*rightValues;
                values = obj.interpolateSpace(values,sourceGrid,isVelocity);
                series.setSnapshot(nextIndex,values);
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
                error('AcousticSolver1DStaggered:SpatialInterpolationFailed', ...
                    'Spatial interpolation produced a nonfinite output value.');
            end
            values = values(:);
        end

        function values = dataOnGrid(~,data,grid,fieldName)
            target = GridFunction(grid);

            if isa(data,'GridFunction')
                if data.dim ~= grid.dim
                    error('AcousticSolver1DStaggered:DimensionMismatch', ...
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
                catch cause
                    error('AcousticSolver1DStaggered:InvalidInitialData', ...
                        'Invalid initial %s data: %s',fieldName,cause.message);
                end
                values = target.values;
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
