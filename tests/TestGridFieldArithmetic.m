classdef TestGridFieldArithmetic < matlab.unittest.TestCase
    methods (Test)
        function timePointsAndTimesOverloadCoexist(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],3));
            U = GridFunctionTimeSeries(grid,[0;0.5;1]);
            U.setValues(@(x,t) x+t);
            V = times(U,U);
            testCase.verifyEqual(U.timePoints,[0,0.5,1]);
            testCase.verifyEqual(U.Nt,3);
            testCase.verifyEqual(V.timePoints,U.timePoints);
            testCase.verifyEqual(V.values,U.values.^2);
        end

        function supportedOperationsAndIndependentResults(testCase)
            x = Grid1DUnifPrimal([0,1],3);
            y = Grid1DUnifPrimal([0,1],2);
            point = Grid1D([0,1],0.5);
            grids = {GridSpace1D(x),GridSpace2D(x,y), ...
                ReceiverSpace1D([0.8;0.2;0.8]), ...
                ReceiverSpace2D([0.8,1;0.2,0]), ...
                GridSpace2D(x,point)};
            for j = 1:numel(grids)
                for series = [false,true]
                    if series
                        a = GridFunctionTimeSeries(grids{j},[0,0.25,1], ...
                            Label="a",Units="Pa",TimeLabel="t",TimeUnits="s");
                        b = GridFunctionTimeSeries(grids{j},[0,0.25,1], ...
                            Label="b",Units="Pa",TimeLabel="t",TimeUnits="s");
                    else
                        a = GridFunction(grids{j},Label="a",Units="Pa");
                        b = GridFunction(grids{j},Label="b",Units="Pa");
                    end
                    av = reshape(1:numel(a.values),size(a.values)) + 1i;
                    bv = 2*ones(size(b.values));
                    a.setValues(av);
                    b.setValues(bv);
                    actual = {a+b,a-b,-a,a.*b,a./b, ...
                        3*a,a*3,3.*a,a.*3,a./3,3./a, ...
                        a+3,3+a,a-3,3-a};
                    expected = {av+bv,av-bv,-av,av.*bv,av./bv, ...
                        3.*av,av.*3,3.*av,av.*3,av./3,3./av, ...
                        av+3,3+av,av-3,3-av};
                    for k = 1:numel(actual)
                        r = actual{k};
                        testCase.verifyClass(r,class(a));
                        testCase.verifyEqual(r.grid,a.grid);
                        testCase.verifyEqual(r.values,expected{k});
                        if series
                            testCase.verifyEqual(r.timePoints,a.timePoints);
                            testCase.verifyEqual(r.timeLabel,"t");
                            testCase.verifyEqual(r.timeUnits,"s");
                        end
                        r.setValues(zeros(size(r.values)));
                        testCase.verifyEqual(a.values,av);
                        testCase.verifyEqual(b.values,bv);
                    end
                    product = a.*b;
                    quotient = a./b;
                    inverse = 3./a;
                    sumField = a+b;
                    scaled = 3*a;
                    testCase.verifyEqual(product.units,"(Pa)*(Pa)");
                    testCase.verifyEqual(quotient.units,"");
                    testCase.verifyEqual(inverse.units,"1/(Pa)");
                    testCase.verifyEqual(sumField.units,"Pa");
                    testCase.verifyEqual(sumField.label,"(a) + (b)");
                    testCase.verifyEqual(scaled.label,"a");
                    testCase.verifyEqual(scaled.units,"Pa");
                end
            end
        end

        function incompatibleOperandsAreRejected(testCase)
            grid = GridSpace1D(Grid1DUnifPrimal([0,1],3));
            otherGrid = GridSpace1D(Grid1DUnifPrimal([0,2],3));
            a = GridFunction(grid,Units="Pa");
            b = GridFunction(otherGrid,Units="Pa");
            testCase.verifyError(@() a+b,'GridFieldArithmetic:GridMismatch');
            b = GridFunction(grid,Units="m");
            testCase.verifyError(@() a+b,'GridFieldArithmetic:UnitsMismatch');
            testCase.verifyError(@() a-b,'GridFieldArithmetic:UnitsMismatch');
            testCase.verifyError(@() a*b, ...
                'GridFieldArithmetic:InvalidMultiplication');
            testCase.verifyError(@() a+[1,2,3], ...
                'GridFieldArithmetic:InvalidOperand');
            testCase.verifyError(@() [1,2,3].*a, ...
                'GridFieldArithmetic:InvalidOperand');
            testCase.verifyError(@() a+"bad", ...
                'GridFieldArithmetic:InvalidOperand');
            U = GridFunctionTimeSeries(grid,[0,1],TimeUnits="s");
            V = GridFunctionTimeSeries(grid,[0,2],TimeUnits="s");
            testCase.verifyError(@() U-V,'GridFieldArithmetic:TimeMismatch');
            V = GridFunctionTimeSeries(grid,[0,1],TimeUnits="ms");
            testCase.verifyError(@() U+V,'GridFieldArithmetic:TimeUnitsMismatch');
            testCase.verifyError(@() a+U,'GridFieldArithmetic:ClassMismatch');
            testCase.verifyError(@() U.*a,'GridFieldArithmetic:ClassMismatch');
        end

        function receiverOrderMustMatch(testCase)
            a = GridFunction(ReceiverSpace1D([0;1]));
            b = GridFunction(ReceiverSpace1D([1;0]));
            testCase.verifyError(@() a-b,'GridFieldArithmetic:GridMismatch');
        end

        function singletonTimeAndUnitLabels(testCase)
            grid = GridSpace2D(Grid1D([0,1],0.5),Grid1D([0,1],0.5));
            a = GridFunctionTimeSeries(grid,0,Units="Pa");
            b = GridFunctionTimeSeries(grid,0,Units="s");
            a.setValues(6);
            b.setValues(2);
            r = a./b;
            testCase.verifyEqual(r.values,3);
            testCase.verifyEqual(r.units,"(Pa)/(s)");
            b.units = "";
            r = a.*b;
            testCase.verifyEqual(r.units,"Pa");
            r = a./b;
            testCase.verifyEqual(r.units,"Pa");
        end
    end
end
