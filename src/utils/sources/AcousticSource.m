classdef (Abstract) AcousticSource < handle
    % AcousticSource is the common base class for acoustic source objects.
    %
    % Concrete subclasses define the spatial dimension and acoustic field
    % components. The base class provides a common type for solver and
    % source-assembly code.

    methods (Abstract)
        append(obj,other)
        % Append all terms from a compatible acoustic source.
    end
end
