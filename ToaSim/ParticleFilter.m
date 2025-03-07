classdef ParticleFilter < BaseFilter
    properties
        Property1
    end

    methods
        function obj = ParticleFilter(inputArg1)
            obj@BaseFilter(inputArg1);
            obj.Property1 = inputArg1^2;
            disp(obj.Property1);
            
        end
    end
end