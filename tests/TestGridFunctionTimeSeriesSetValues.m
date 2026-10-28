classdef TestGridFunctionTimeSeriesSetValues < matlab.unittest.TestCase
    methods (Test)
        function arraysAndFunctionsOnEverySpace(testCase)
            x = Grid1DUnifPrimal([0,1],5);
            y = Grid1DUnifPrimal([0,2],3);
            point = Grid1D([0,1],0.5);
            grids = {GridSpace1D(x), GridSpace2D(x,y), ...
                ReceiverSpace1D([0.8;0.2;0.8]), ...
                ReceiverSpace2D([0.8,1;0.2,0;0.8,1]), ...
                GridSpace1D(point), GridSpace2D(point,y), ...
                GridSpace2D(x,point), GridSpace2D(point,point)};
            timeGrids = {[0,0.25,1], 0.25};

            for j = 1:numel(grids)
                grid = grids{j};
                for i = 1:numel(timeGrids)
                    times = timeGrids{i};
                    U = GridFunctionTimeSeries(grid,times, ...
                        Label="Pressure",Units="Pa",TimeUnits="s");
                    if grid.dim == 1
                        f = @(x,t) (1+2i)*x + 3*t;
                        spatial = @(x) (1+2i)*x;
                    else
                        f = @(x,y,t) (1+2i)*x - 2*y + 3*t;
                        spatial = @(x,y) (1+2i)*x - 2*y;
                    end
                    U.setValues(f);
                    base = GridFunction(grid);
                    base.setValues(spatial);
                    for k = 1:U.Nt
                        u = U.snapshot(k);
                        testCase.verifyEqual(u.values, ...
                            base.values+3*times(k),AbsTol=1e-14);
                    end
                    expected = U.values;
                    U.setValues(zeros(size(expected)));
                    U.setValues(expected);
                    testCase.verifyEqual(U.values,expected);
                    testCase.verifyEqual(U.times,times);
                    testCase.verifyEqual(U.label,"Pressure");
                    testCase.verifyEqual(U.units,"Pa");
                    testCase.verifyEqual(U.timeUnits,"s");
                end
            end
        end

        function rejectsInvalidInputsWithoutChangingValues(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],5));
            U = GridFunctionTimeSeries(grid,[0,0.25,1]);
            expected = reshape(1:15,5,3);
            U.setValues(expected);
            testCase.verifyError(@() U.setValues(expected.'), ...
                'GridFunctionTimeSeries:InvalidSize');
            testCase.verifyError(@() U.setValues(expected(:)), ...
                'GridFunctionTimeSeries:InvalidSize');
            testCase.verifyError(@() U.setValues("invalid"), ...
                'GridFunctionTimeSeries:InvalidData');
            testCase.verifyError(@() U.setValues(@(x,t) "invalid"), ...
                'GridFunctionTimeSeries:InvalidData');
            testCase.verifyError(@() U.setValues(@(x,t) x(1:end-1)), ...
                'GridFunctionTimeSeries:InvalidSize');
            testCase.verifyEqual(U.values,expected);
            testCase.verifyError(@() U.setValues(@failAtLaterTime), ...
                'TestGridFunctionTimeSeriesSetValues:EvaluationFailed');
            testCase.verifyEqual(U.values,expected);
        end

        function snapshotStillUpdatesOnlyOneTime(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],5));
            U = GridFunctionTimeSeries(grid,[0,0.25,1]);
            U.setValues(@(x,t) x+t);
            expected = U.values;
            expected(:,2) = 9;
            U.setSnapshot(2,9*ones(5,1));
            testCase.verifyEqual(U.values,expected);
        end
    end
end

function values = failAtLaterTime(x,t)
    if t > 0
        error('TestGridFunctionTimeSeriesSetValues:EvaluationFailed', ...
            'Deliberate failure after the first snapshot.');
    end
    values = x;
end
