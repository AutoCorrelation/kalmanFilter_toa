classdef Environment
    properties
        Property1
    end

    methods
        function obj = Environment(inputArg1)
            obj.Property1 = inputArg1^2;
            disp(obj.Property1);
            
        end
    end
end