classdef TestGridFieldMetadata < matlab.unittest.TestCase
    methods (Test)
        function copiesSharedMetadataAcrossBothClasses(testCase)
            % Different grids/dimensions are allowed for metadata transfer.
            sourceGrid = GridSpace1D(Grid1DUnifPrimal([0,1],3));
            targetGrid = ReceiverSpace2D([0.2,0.4;0.8,0.1]);
            for sourceIsSeries = [false,true]
                for targetIsSeries = [false,true]
                    if sourceIsSeries
                        source = GridFunctionTimeSeries(sourceGrid,[0,1], ...
                            TimeLabel="Source time",TimeUnits="ms");
                    else
                        source = GridFunction(sourceGrid);
                    end
                    source.label = "Pressure";
                    source.units = "Pa";
                    sourceValues = 7*ones(size(source.values));
                    source.setValues(sourceValues);
                    if targetIsSeries
                        target = GridFunctionTimeSeries(targetGrid,[2,3,4], ...
                            TimeLabel="Target time",TimeUnits="s");
                    else
                        target = GridFunction(targetGrid);
                    end
                    targetValues = reshape(1:numel(target.values),size(target.values));
                    target.setValues(targetValues);
                    target.copyMetadataFrom(source);
                    testCase.verifyEqual(target.label,"Pressure");
                    testCase.verifyEqual(target.units,"Pa");
                    testCase.verifyEqual(target.values,targetValues);
                    testCase.verifyEqual(target.grid,targetGrid);
                    testCase.verifyEqual(source.values,sourceValues);
                    if targetIsSeries
                        testCase.verifyEqual(target.timePoints,[2,3,4]);
                        if sourceIsSeries
                            testCase.verifyEqual(target.timeLabel,"Source time");
                            testCase.verifyEqual(target.timeUnits,"ms");
                        else
                            testCase.verifyEqual(target.timeLabel,"Target time");
                            testCase.verifyEqual(target.timeUnits,"s");
                        end
                    end
                    % Metadata strings are copied, not linked to the source.
                    source.label = "Changed source";
                    source.units = "kPa";
                    testCase.verifyEqual(target.label,"Pressure");
                    testCase.verifyEqual(target.units,"Pa");
                end
            end
        end

        function emptyMetadataAndSelfCopy(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],3));
            source = GridFunctionTimeSeries(grid,[0,1], ...
                Label="",Units="",TimeLabel="",TimeUnits="");
            target = GridFunctionTimeSeries(grid,[2,3]);
            target.setValues(5*ones(size(target.values)));
            target.copyMetadataFrom(source);
            target.copyMetadataFrom(target);
            testCase.verifyEqual(target.label,"");
            testCase.verifyEqual(target.units,"");
            testCase.verifyEqual(target.timeLabel,"");
            testCase.verifyEqual(target.timeUnits,"");
            testCase.verifyEqual(target.values,5*ones(3,2));
            testCase.verifyEqual(target.timePoints,[2,3]);
        end

        function invalidSourcesLeaveMetadataUnchanged(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],3));
            targets = {GridFunction(grid,Label="Keep",Units="Pa"), ...
                GridFunctionTimeSeries(grid,[0,1],Label="Keep",Units="Pa")};
            invalid = {[],42,"Pressure",struct('label',"x",'units',"m")};
            for j = 1:numel(targets)
                target = targets{j};
                for k = 1:numel(invalid)
                    testCase.verifyError(@() target.copyMetadataFrom(invalid{k}), ...
                        [class(target),':InvalidMetadataSource']);
                    testCase.verifyEqual(target.label,"Keep");
                    testCase.verifyEqual(target.units,"Pa");
                end
            end
        end
    end
end
