classdef TestGridFieldNorm < matlab.unittest.TestCase
    methods (Test)
        function vectorNormsAndAllReductionModes(testCase)
            grid = GridSpace2D(Grid1DUnifPrimal([0,2],3), ...
                Grid1DUnifPrimal([0,1],2));
            U = GridFunctionTimeSeries(grid,[0,0.2,1], ...
                Label="Pressure",Units="Pa",TimeUnits="s");
            v = reshape(1:18,3,2,3) + 1i;
            U.setValues(v);
            original = U.values;
            for p = [1,2,Inf]
                testCase.verifyEqual(norm(U,p),norm(v(:),p),RelTol=1e-13);
                spatial = norm(U,p,Over="space");
                temporal = norm(U,p,Over="time");
                firstAxis = norm(U,p,Dimension=1);
                mixed = norm(U,p,Dimension=[1,3]);
                testCase.verifyClass(spatial,'GridFunctionTimeSeries');
                testCase.verifyEqual(spatial.N,[1,1]);
                testCase.verifyEqual(spatial.timePoints,U.timePoints);
                testCase.verifyEqual(spatial.timeUnits,"s");
                testCase.verifyClass(temporal,'GridFunction');
                testCase.verifyEqual(temporal.grid,U.grid);
                testCase.verifyEqual(temporal.units,"Pa");
                testCase.verifyEqual(firstAxis.N,[1,2]);
                testCase.verifyEqual(firstAxis.grid.y,grid.y);
                testCase.verifyEqual(firstAxis.grid.x.pts,1);
                testCase.verifyEqual(mixed.N,[1,2]);
                for k = 1:3
                    snapshot = v(:,:,k);
                    testCase.verifyEqual(spatial.values(1,1,k), ...
                        norm(snapshot(:),p),RelTol=1e-13);
                    for j = 1:2
                        testCase.verifyEqual(firstAxis.values(1,j,k), ...
                            norm(v(:,j,k),p),RelTol=1e-13);
                    end
                end
                for j = 1:2
                    slice = v(:,j,:);
                    testCase.verifyEqual(mixed.values(1,j), ...
                        norm(slice(:),p),RelTol=1e-13);
                    for i = 1:3
                        testCase.verifyEqual(temporal.values(i,j), ...
                            norm(reshape(v(i,j,:),[],1),p),RelTol=1e-13);
                    end
                end
            end
            testCase.verifyEqual(U.values,original);
        end

        function automaticIntegralQuadrature(testCase)
            grid = GridSpace2D(Grid1DUnifPrimal([0,2],5), ...
                Grid1DUnifPrimal([0,3],4));
            U = GridFunctionTimeSeries(grid,[0,0.25,1], ...
                Units="Pa",TimeUnits="s");
            U.setValues(2*ones(size(U.values)));
            testCase.verifyEqual(norm(U,2,Weighted=true),2*sqrt(6),AbsTol=1e-13);
            spatial = norm(U,2,Over="space",Weighted=true,MeasureUnits="m^2");
            testCase.verifyEqual(spatial.values, ...
                2*sqrt(6)*ones(1,1,3),AbsTol=1e-13);
            testCase.verifyEqual(spatial.units,"(Pa)*(m^2)^(1/2)");
            temporal = norm(U,2,Over="time",Weighted=true,MeasureUnits="s");
            testCase.verifyEqual(temporal.values,2*ones(5,4),AbsTol=1e-13);
            testCase.verifyEqual(temporal.units,"(Pa)*(s)^(1/2)");
            U.setValues(@(x,y,t) (1+t)*ones(size(x)));
            temporal = norm(U,2,Over="time",Weighted=true);
            expected = sqrt(trapz(U.timePoints,(1+U.timePoints).^2));
            testCase.verifyEqual(temporal.values,expected*ones(5,4),AbsTol=1e-13);
            dual = GridSpace1D(Grid1DUnifDual(grid.x));
            u = GridFunction(dual);
            u.setValues(3*ones(size(u.values)));
            testCase.verifyEqual(norm(u,2,Weighted=true),3*sqrt(2),AbsTol=1e-13);
            % A one-cell dual grid still has a known cell width.
            oneCell = GridSpace1D(Grid1DUnifDual(Grid1DUnifPrimal([0,2],2)));
            u = GridFunction(oneCell);
            u.setValues(3);
            testCase.verifyEqual(norm(u,2,Weighted=true),3*sqrt(2),AbsTol=1e-13);
        end

        function explicitTensorAndSeparableWeights(testCase)
            grid = GridSpace2D(Grid1DUnifPrimal([0,1],3), ...
                Grid1DUnifPrimal([0,1],2));
            U = GridFunctionTimeSeries(grid,[0,0.5,1]);
            v = reshape(1:18,3,2,3);
            U.setValues(v);
            wx = [1;2;1];
            wy = [2;3];
            wt = [1;2;4];
            w = reshape(wx,3,1,1).*reshape(wy,1,2,1).*reshape(wt,1,1,3);
            expected = sqrt(sum(w(:).*v(:).^2));
            testCase.verifyEqual(norm(U,2,Weights=w),expected,RelTol=1e-13);
            testCase.verifyEqual(norm(U,2,Weights={wx,wy,wt}),expected,RelTol=1e-13);
            r = norm(U,2,Dimension=[1,3],Weights={wx,wt});
            for j = 1:2
                slice = reshape(v(:,j,:),3,3);
                weights = wx*wt.';
                testCase.verifyEqual(r.values(1,j), ...
                    sqrt(sum(weights(:).*slice(:).^2)),RelTol=1e-13);
            end
        end

        function receiverWeightsAndInfinitySupport(testCase)
            grid = ReceiverSpace2D([0,1;1,0;2,2]);
            U = GridFunctionTimeSeries(grid,[0,1]);
            U.setValues([3,4;100,200;4,3]);
            w = [2;0;1];
            r = norm(U,2,Over="space",Weights=w);
            testCase.verifyEqual(r.values,sqrt([34,41]),AbsTol=1e-13);
            testCase.verifyEqual(r.grid.points,mean(grid.points,1));
            r = norm(U,Inf,Over="space",Weights=w);
            testCase.verifyEqual(r.values,[4,4]);
            r = norm(U,2,Over="time",Weighted=true);
            testCase.verifyEqual(r.values,sqrt([12.5;25000;12.5]),AbsTol=1e-12);
            testCase.verifyError(@() norm(U,2,Over="space",Weighted=true), ...
                'GridFieldNorm:ReceiverWeightsRequired');
        end

        function singletonAxesZeroValuesAndStableScaling(testCase)
            grid = GridSpace1D(Grid1D([0,1],0.5));
            U = GridFunctionTimeSeries(grid,2);
            U.setValues(3);
            testCase.verifyEqual(norm(U),3);
            r = norm(U,2,Over="space");
            testCase.verifyEqual(r.values,3);
            r = norm(U,2,Over="time");
            testCase.verifyEqual(r.values,3);
            testCase.verifyError(@() norm(U,2,Weighted=true), ...
                'GridFieldNorm:SingletonWeightsRequired');
            testCase.verifyEqual(norm(U,2,Weights={2,4}),sqrt(72),AbsTol=1e-13);
            U.setValues(0);
            testCase.verifyEqual(norm(U),0);
            u = GridFunction(GridSpace1D(Grid1DUnifPrimal([0,1],2)));
            u.setValues([1e200;1e200]);
            testCase.verifyEqual(norm(u),sqrt(2)*1e200,RelTol=1e-13);
            u.setValues([NaN;2]);
            testCase.verifyTrue(isnan(norm(u)));
            testCase.verifyEqual(norm(u,2,Weights=[0;1]),2);
        end

        function spatialFieldAndLineReductions(testCase)
            x = Grid1DUnifPrimal([0,1],2);
            u = GridFunction(GridSpace2D(x,x));
            u.setValues(eye(2));
            testCase.verifyEqual(norm(u,2),sqrt(2),AbsTol=1e-14);
            r = norm(u,2,Dimension=1);
            testCase.verifyClass(r,'GridFunction');
            testCase.verifyEqual(r.values,[1,1]);
            line = GridSpace2D(Grid1D([0,1],0.25), ...
                Grid1DUnifPrimal([0,2],3));
            u = GridFunction(line);
            u.setValues(3*ones(1,3));
            r = norm(u,2,Dimension=2,Weighted=true);
            testCase.verifyEqual(r.values,3*sqrt(2),AbsTol=1e-13);
            testCase.verifyEqual(r.grid.x,line.x);
            % Scaling must occur after incorporating weight magnitudes.
            u = GridFunction(GridSpace1D(x));
            u.setValues([1e200;1]);
            testCase.verifyEqual(norm(u,2,Weights=[1e-300;1e300]), ...
                1e150,RelTol=1e-13);
        end

        function rejectsInvalidOptions(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],3));
            u = GridFunction(grid);
            U = GridFunctionTimeSeries(grid,[0,1]);
            testCase.verifyError(@() norm(u,3),'GridFieldNorm:InvalidOrder');
            testCase.verifyError(@() norm(u,2,Over="time"),'GridFieldNorm:NoTimeAxis');
            testCase.verifyError(@() norm(U,2,Over="bad"),'GridFieldNorm:InvalidReduction');
            testCase.verifyError(@() norm(U,2,Over="space",Dimension=1), ...
                'GridFieldNorm:ConflictingReduction');
            for dims = {0,3,[1,1],1.5,NaN}
                testCase.verifyError(@() norm(U,2,Dimension=dims{1}), ...
                    'GridFieldNorm:InvalidDimension');
            end
            for w = {[-1;1;1],zeros(3,1),[NaN;1;1],ones(4,1),"bad"}
                testCase.verifyError(@() norm(U,2,Over="space",Weights=w{1}), ...
                    'GridFieldNorm:InvalidWeights');
            end
        end
    end
end
