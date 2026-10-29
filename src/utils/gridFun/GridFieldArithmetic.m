classdef GridFieldArithmetic
    % Shared implementation of GridFunction and time-series arithmetic.
    % Numeric scalars are interpreted in the field's units for + and -,
    % and as dimensionless factors for multiplication and division.
    % Unit expressions are labels only; no conversion is performed.
    % Mixed GridFunction/time-series arithmetic is intentionally unsupported.

    methods (Static)
        function result = binary(a,b,operation)
            fieldA = isa(a,'GridFunction') || isa(a,'GridFunctionTimeSeries');
            fieldB = isa(b,'GridFunction') || isa(b,'GridFunctionTimeSeries');
            validA = isscalar(a) && (fieldA || isnumeric(a));
            validB = isscalar(b) && (fieldB || isnumeric(b));
            if ~validA || ~validB || ~(fieldA || fieldB)
                error('GridFieldArithmetic:InvalidOperand', ...
                    'Operands must be scalar field objects or numeric scalars.');
            end

            if fieldA
                source = a;
                valuesA = a.values;
            else
                source = b;
                valuesA = a;
            end
            if fieldB
                valuesB = b.values;
            else
                valuesB = b;
            end

            if operation == "mtimes" && fieldA && fieldB
                error('GridFieldArithmetic:InvalidMultiplication', ...
                    'Use .* for pointwise multiplication of two fields.');
            end
            if fieldA && fieldB
                if ~strcmp(class(a),class(b))
                    error('GridFieldArithmetic:ClassMismatch', ...
                        'Both field operands must have the same class.');
                end
                if ~isequal(a.grid,b.grid)
                    error('GridFieldArithmetic:GridMismatch', ...
                        'Both field operands must have matching spatial grids.');
                end
                if isa(a,'GridFunctionTimeSeries')
                    if ~isequal(a.timePoints,b.timePoints)
                        error('GridFieldArithmetic:TimeMismatch', ...
                            'Both time series must have identical stored times.');
                    end
                    if a.timeUnits ~= b.timeUnits
                        error('GridFieldArithmetic:TimeUnitsMismatch', ...
                            'Both time series must have identical time-unit strings.');
                    end
                end
            end

            units = source.units;
            if operation == "plus" || operation == "minus"
                if fieldA && fieldB && a.units ~= b.units
                    error('GridFieldArithmetic:UnitsMismatch', ...
                        'Addition and subtraction require identical unit strings.');
                end
            elseif operation == "times" || operation == "rdivide"
                unitsA = "";
                unitsB = "";
                if fieldA
                    unitsA = a.units;
                end
                if fieldB
                    unitsB = b.units;
                end
                units = GridFieldArithmetic.combineUnits( ...
                    unitsA,unitsB,operation);
            end

            switch operation
                case "plus"
                    values = valuesA + valuesB;
                    symbol = " + ";
                case "minus"
                    values = valuesA - valuesB;
                    symbol = " - ";
                case {"times","mtimes"}
                    values = valuesA .* valuesB;
                    symbol = " .* ";
                case "rdivide"
                    values = valuesA ./ valuesB;
                    symbol = " ./ ";
                otherwise
                    error('GridFieldArithmetic:InvalidOperation', ...
                        'Unsupported arithmetic operation.');
            end

            label = source.label;
            if fieldA && fieldB
                label = "(" + a.label + ")" + symbol + "(" + b.label + ")";
            end
            if isa(source,'GridFunctionTimeSeries')
                result = GridFunctionTimeSeries(source.grid,source.timePoints, ...
                    Label=label,Units=units, ...
                    TimeLabel=source.timeLabel,TimeUnits=source.timeUnits);
            else
                result = GridFunction(source.grid,Label=label,Units=units);
            end
            result.setValues(values);
        end
    end

    methods (Static, Access = private)
        function units = combineUnits(a,b,operation)
            if operation == "times"
                if a == ""
                    units = b;
                elseif b == ""
                    units = a;
                else
                    units = "(" + a + ")*(" + b + ")";
                end
            else
                if b == ""
                    units = a;
                elseif a == b
                    units = "";
                elseif a == ""
                    units = "1/(" + b + ")";
                else
                    units = "(" + a + ")/(" + b + ")";
                end
            end
        end
    end
end
