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
        P0
        R
        kf1x
        kf1P
    end

    methods
        function obj = simulate(iteration, numPoints)
            clc;

            obj.iteration = iteration;
            obj.numPoints = numPoints;
            obj.noiseVariance = [0.01; 0.1; 1; 10; 100];
            obj.A = eye(2);
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
            toaError = pos - obj.toaPos(:, :, step, :);
            obj.w_bias = squeeze(mean(errorPos, 2));    % 2x5
            Eerror = squeeze(mean(toaError, 2));
            x = squeeze(errorPos);
            toaError = squeeze(toaError);
            eeT = zeros(2, 2, obj.iteration, length(obj.noiseVariance));
            xxT = zeros(2, 2, obj.iteration, length(obj.noiseVariance));
            for iter = 1:obj.iteration
                for noise = 1:length(obj.noiseVariance)
                    eeT(:,:,iter,noise) = x(:,iter,noise)*x(:,iter,noise)';

                    xxT(:,:,iter,noise) = toaError(:,iter,noise)*toaError(:,iter,noise)';
                end
            end
            EeeT = squeeze(mean(eeT,3));
            ExxT = squeeze(mean(xxT,3));
            for noise = 1:length(obj.noiseVariance)
                obj.Q(:, :, noise) = EeeT(:,:,noise) - obj.w_bias(:, noise) * obj.w_bias(:, noise)';

                obj.P0(:, :, noise) = ExxT(:, :, noise) - Eerror(:, noise) * Eerror(:, noise)';
            end
        end

        function obj = getR(obj)
            obj.R = zeros(6, 6, obj.numPoints, length(obj.noiseVariance));
            for noise = 1:length(obj.noiseVariance)
                for step = 1:obj.numPoints
                    for iter = 1:obj.iteration
                        obj.R(:, :, step, noise) = obj.R(:, :, step, noise) + ...
                            [4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 0;...
                            4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 0 -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2);...
                            4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(1,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2);...
                            -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2+obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2);...
                            -4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 0 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(2,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2);...
                            0 -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) -4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(4,iter,step,noise)^2) 4*obj.noiseVariance(noise)*(obj.rangingInfo(3,iter,step,noise)^2+obj.rangingInfo(4,iter,step,noise)^2)
                            ];
                    end
                    obj.R(:, :, step, noise) = obj.R(:, :, step, noise) / obj.iteration;
                end
            end
        end

        function obj = ToaKalmanFilter(obj)
            obj.kf1x = zeros(2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            obj.kf1P = zeros(2, 2, obj.iteration, obj.numPoints, length(obj.noiseVariance));
            for iter = 1:obj.iteration
                for step = 1:obj.numPoints
                    for noise = 1:length(obj.noiseVariance)
                        switch step
                            case 1
                                obj.kf1x(:,iter,step,noise) = [0; 0];
                                obj.kf1P(:,:,iter,step,noise) = obj.P0(:, :, noise);
                                velocity = 0;
                            case 2
                                obj.kf1x(:,iter,step,noise) = obj.toaPos(:, iter, step, noise);
                                obj.kf1P(:,:,iter,step,noise) = obj.kf1P(:,:,iter,step-1,noise);
                            case 3
                                obj.kf1x(:,iter,step,noise) = obj.toaPos(:, iter, step, noise);
                                obj.kf1P(:,:,iter,step,noise) = obj.kf1P(:,:,iter,step-1,noise);
                                velocity = (obj.kf1x(:,iter,step,noise) - obj.kf1x(:,iter,step-1,noise)) / 0.1;
                            otherwise
                                [obj.kf1x(:,iter,step,noise), obj.kf1P(:,:,iter,step,noise)] = ToaKf(obj.kf1x(:,iter,step-1,noise), obj.kf1P(:,:,iter,step-1,noise), 0, velocity, obj.A, obj.Q(:, :, noise), obj.w_bias(:, noise), obj.H, obj.R(:, :, step, noise), obj.z(:, iter, step, noise));
                                velocity = (obj.kf1x(:,iter,step,noise) - obj.kf1x(:,iter,step-1,noise)) / 0.1;
                        end
                    end
                end
            end
        end

        function obj = resultPlot(obj)
            toaRMSE = zeros(length(obj.noiseVariance),1);
            kf1RMSE = zeros(length(obj.noiseVariance),1);
            for iter = 1:obj.iteration
                for step = 2:obj.numPoints
                    pos = [step-1; step-1];
                    for noise = 1:length(obj.noiseVariance)
                        toaRMSE(noise,1) = toaRMSE(noise,1) + norm(pos - obj.toaPos(:, iter, step, noise));
                        kf1RMSE(noise,1) = kf1RMSE(noise,1) + norm(pos - obj.kf1x(:, iter, step, noise));

                    end
                end
            end

            toaRMSE = toaRMSE / (obj.iteration * (obj.numPoints-1));
            kf1RMSE = kf1RMSE / (obj.iteration * (obj.numPoints-1));


            figure;
            semilogx(obj.noiseVariance, toaRMSE, 'o-','DisplayName','TOA');
            hold on;
            semilogx(obj.noiseVariance, kf1RMSE, 'o-','DisplayName','KF1');
            xlabel('Noise Variance');
            ylabel('RMSE');
            legend;
        end
    end
end