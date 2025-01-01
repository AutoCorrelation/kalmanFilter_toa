classdef simulate
    properties
        iteration
        numPoints
        anchorPos
        rangingInfo
        noiseVariance
    end

    methods
        function obj = simulate(iteration, numPoints)
            clc;

            obj.iteration = iteration;
            obj.numPoints = numPoints;
            obj.noiseVariance = [0.01; 0.1; 1; 10; 100];
        end

        function obj = setsquareAnchorPos(obj, d)
            obj.anchorPos = [0, d; 0, 0; d, 0; d, d];
        end

        function obj = generateRangingInfo(obj)
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        exactPos = [step-1; step-1];
                        for achor = 1:4
                            obj.rangingInfo(iter, step, noise, achor) ...
                                = norm(exactPos - obj.anchorPos(achor, :)) ...
                                + sqrt(obj.noiseVariance(noise)) * randn;
                        end
                    end
                end
            end
        end

        
    end
end