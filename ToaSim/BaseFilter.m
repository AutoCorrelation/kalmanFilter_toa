classdef (Abstract) BaseFilter
    properties
        Property1
    end

    methods (Abstract)
        function obj = BaseFilter(inputArg1)
            obj.Property1 = inputArg1^2;
            disp(obj.Property1);
            
        end
    end
end