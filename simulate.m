classdef simulate
    properties
        iteration
        numPoints
        anchorPos
        rangingInfo
        noiseVariance
        d
        z
        H
        toaPos
        A
        Q
        w_bias
    end

    methods
        function obj = simulate(iteration, numPoints)
            clc;

            obj.iteration = iteration;
            obj.numPoints = numPoints;
            obj.noiseVariance = [0.01; 0.1; 1; 10; 100];
        end

        function obj = setsquareAnchorPos(obj, d)
            obj.d = d;
            obj.anchorPos = [0, d; 0, 0; d, 0; d, d];
        end

        function obj = generateRangingInfo(obj)
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        exactPos = [step-1; step-1];
                        for anchor = 1:4
                            obj.rangingInfo(anchor, iter, step, noise) ...
                                = norm(exactPos - obj.anchorPos(anchor, :)') ...
                                + sqrt(obj.noiseVariance(noise)) * randn;
                        end
                    end
                end
            end
        end

        function obj = toa(obj)
            obj.H = [...
                0, -2*obj.d
                2*obj.d, -2*obj.d
                2*obj.d, 0
                2*obj.d, 0
                2*obj.d, 2*obj.d
                0, 2*obj.d];
            pseudoInverseH = (obj.H' * obj.H) \ obj.H';
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        zvec = [obj.rangingInfo(1, iter, step, noise)^2 - obj.rangingInfo(2, iter, step, noise)^2 - obj.d^2;
                            obj.rangingInfo(1, iter, step, noise)^2 - obj.rangingInfo(3, iter, step, noise)^2;
                            obj.rangingInfo(1, iter, step, noise)^2 - obj.rangingInfo(4, iter, step, noise)^2 + obj.d^2;
                            obj.rangingInfo(2, iter, step, noise)^2 - obj.rangingInfo(3, iter, step, noise)^2 + obj.d^2;
                            obj.rangingInfo(2, iter, step, noise)^2 - obj.rangingInfo(4, iter, step, noise)^2 + 2*(obj.d^2);
                            obj.rangingInfo(3, iter, step, noise)^2 - obj.rangingInfo(4, iter, step, noise)^2 + obj.d^2];
                        obj.z(:, iter, step, noise) = zvec;
                        obj.toaPos(:, iter, step, noise) = pseudoInverseH * obj.z(:, iter, step, noise);
                    end
                end
            end
        end

        function obj = getQ(obj)
            obj.Q = zeros(2, 2, length(obj.noiseVariance));
            step = 4;
            pos = step * ones(size(obj.toaPos(:, :, step, :)));
            toaVel = (obj.toaPos(:, :, step-1, :) - obj.toaPos(:, :, step-2, :)) / 0.1;   % pos 2, 1
            errorPos = pos - obj.toaPos(:, :, step, :) - toaVel * 0.1;              % pos 
            obj.w_bias = squeeze(mean(errorPos, 2));    % 2x5
            x = squeeze(errorPos);
            for iter = 1:obj.iteration
                for noise = 1:length(obj.noiseVariance)
                    xxT(:,:,iter,noise) = x(:,iter,noise)*x(:,iter,noise)';
                end
            end
            ExxT = squeeze(mean(xxT,3));
            for noise = 1:length(obj.noiseVariance)
                obj.Q(:, :, noise) = ExxT(:,:,noise) - obj.w_bias(:, noise) * obj.w_bias(:, noise)';
            end
            
        end


    end
end